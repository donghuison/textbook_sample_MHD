! Directional split : unsplit, split
! Time marching     : Euler, Predictor-corrector (PC2), Runge-Kutta (RK2), Runge-Kutta (RK3)
! Reconstruction    : NON, MUSCL2, MUSCL3
! Limiter           : MINMOD, SuperBee
#include "config.h"
module timestep
   use parameter
   implicit none
   private
   public :: step_full
contains
   ! ---------------------------------------------------------------------------
   !
   ! ---------------------------------------------------------------------------
   subroutine step_full

#if defined(DIRECTIONAL_SPLIT)
      call step_split
#endif

   end subroutine step_full


   ! ---------------------------------------------------------------------------
   ! Macro for step_split_1d
   ! ---------------------------------------------------------------------------
#if defined(TIME_MARCHING_RK2)
#define STEP_SPLIT_1D(n) step_split_rk2_1d(n)
#else
   ERROR
#endif
   ! ---------------------------------------------------------------------------
   ! fractional time step
   ! ---------------------------------------------------------------------------
   subroutine step_split
#ifdef FLUX_SCHEME_MHD
      use grid, only: V, W, F, Dtime
      use flux_eos, only : source_b, v2u, u2v
#endif !FLUX_SCHEME_MHD
      if (NDIM == 1) then
         call STEP_SPLIT_1D(MX)
      elseif (NDIM == 2) then
         call step_split_2d
      else
         print *, '*** error in step_fractional'
      endif
#ifdef FLUX_SCHEME_MHD
      call v2u(V, W)
      call source_b(F, W, Dtime)
      call u2v(W, V)
#endif !FLUX_SCHEME_MHD
   end subroutine step_split

   ! ---------------------------------------------------------------------------
   ! fractional time step for 2D
   ! ---------------------------------------------------------------------------
   subroutine step_split_2d
      integer :: n
      do n = MX, MY
         call STEP_SPLIT_1D( n )
      end do
   end subroutine step_split_2d
   ! ---------------------------------------------------------------------------
   ! Runge-Kutta 2 in one dimension
   ! ---------------------------------------------------------------------------
   subroutine step_split_rk2_1d(ndir)
      use grid
      use flux_eos
      use boundary
      integer,intent(IN) :: ndir
      call v2u(V, U)
      call boundary_fix(V)
      call get_flux_ndir(ndir)
      W = U
      call w_update_ndir(Dtime, ndir)        ! W := U*
      call u2v(W, V)
      call boundary_fix(V)
      call get_flux_ndir(ndir)
      W = (U + W)*0.5d0
      call w_update_ndir(Dtime*0.5d0, ndir)
      call u2v(W, V)
      call boundary_fix(V)
   end subroutine step_split_rk2_1d


   ! ---------------------------------------------------------------------------
   ! update u by flux for all directions
   ! ---------------------------------------------------------------------------
   subroutine w_update( dt )
      use grid
#ifdef FLUX_SCHEME_MHD
      use flux_eos, only : source_b
#endif !FLUX_SCHEME_MHD
      real(kind=DBL_KIND),intent(IN) :: dt
      integer :: n
      do n = MX, MX+NDIM-1
         call w_update_ndir(dt, n)
      end do
#ifdef FLUX_SCHEME_MHD
      call source_b(F, W, dt)
#endif !FLUX_SCHEME_MHD
   end subroutine w_update
   ! ---------------------------------------------------------------------------
   ! update v by flux for each direction
   ! ---------------------------------------------------------------------------
#define SHIFTR( A, NDIM ) cshift((A),  1, (NDIM) + DIMOFFSET)
#define SHIFTL( A, NDIM ) cshift((A), -1, (NDIM) + DIMOFFSET)
   subroutine w_update_ndir (dt, ndir)
      use grid
      real(kind=DBL_KIND),intent(IN) :: dt
      integer,intent(IN) :: ndir
      real(kind=DBL_KIND),dimension(MX:MZ) :: ds
      integer,parameter :: DIMOFFSET = 1-MX
      ds = get_ds()
      W = W - dt*ds(ndir)*(F(:,:,:,:,ndir) - SHIFTL(F(:,:,:,:,ndir),ndir))
   end subroutine w_update_ndir
   ! ---------------------------------------------------------------------------
   ! flux at cell interface for all directions
   ! ---------------------------------------------------------------------------
   subroutine get_flux(bool_muscl)
      integer :: n
      logical,optional :: bool_muscl
      do n = MX, MX+NDIM-1
         call get_flux_ndir(n, bool_muscl=bool_muscl)
      end do
   end subroutine get_flux
   ! ---------------------------------------------------------------------------
   ! flux at cell interface for each direction
   ! ---------------------------------------------------------------------------
#define MINMOD(x, y) (max(0.d0,min((y)*sign(1.d0,(x)),abs(x)))*sign(1.d0,(x)))
#define SUPERBEE(x, y) (sign(1.d0,(y))*max(0.d0, min(sign(1.d0,(y))*BW*(x),abs(y)), min(sign(1.d0,(y))*(x),BW*abs(y))))
#ifdef MUSCL2_LIMITER_MINMOD
#define FLMT(x, y) MINMOD(x, y)
#endif
#ifdef MUSCL2_LIMITER_SUPERBEE
#define FLMT(x, y) SUPERBEE(x, y)
#endif
#ifdef MUSCL2_WO_LIMITER
#define FLMT(x, y) (y)
#endif
   subroutine get_flux_ndir (ndir, bool_muscl)
      use util
      use grid
      use flux_eos
      integer,intent(IN) :: ndir
      logical,optional :: bool_muscl
      real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX) :: f1d, vl, vr
      integer,dimension(MMIN:MMAX) :: mcycle
      integer :: io,jo,ko,i2,j2,k2,i,j,k,m
#ifdef RECONSTRUCTION_MUSCL2
      real(kind=DBL_KIND),parameter :: BW = 2.d0
#endif
      logical :: b_muscl
      b_muscl = .TRUE.
      if (present(bool_muscl)) b_muscl = bool_muscl

      call util_arroffset(ndir,io,jo,ko)
      i2 = io*2
      j2 = jo*2
      k2 = ko*2

      do m = MMIN, MMAX
         do k = KMIN-ko, KMAX
            do j = JMIN-jo, JMAX
               do i = IMIN-io, IMAX
                  vl(i,j,k,m) = V(i,j,k,m)
                  vr(i,j,k,m) = V(i+io,j+jo,k+ko,m)
               end do
            end do
         end do
      end do

#if !defined(RECONSTRUCTION_NONE)
      if (b_muscl) then
         do m = MMIN, MMAX
            do k = KMIN-ko, KMAX
#if defined(RECONSTRUCTION_MUSCL3)
               !$omp parallel do private(i,dva,dvb)
#elif defined(RECONSTRUCTION_LIMO3)
               !$omp parallel do private(i,duLL,duLR,duRL,duRR,thtL,thtR,etaL,etaR,flagL,flagR,phiL,phiR)
#else
               !$omp parallel do private(i)
#endif
               do j = JMIN-jo, JMAX
                  do i = IMIN-io, IMAX
#if defined(RECONSTRUCTION_MUSCL2)
                     vl(i,j,k,m) = vl(i,j,k,m) &
                        + (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i,j,k,m)-V(i-io,j-jo,k-ko,m)))*0.5d0
                     vr(i,j,k,m) = vr(i,j,k,m) &
                        - (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i+i2,j+j2,k+k2,m)-V(i+io,j+jo,k+ko,m)))*0.5d0
#else
                     ERROR
#endif
                  enddo
               enddo
               !$omp end parallel do
            enddo
         enddo
      endif
#endif !RECONSTRUCTION_NONE
      mcycle = cyclecomp( ndir )
      vl = vl(:,:,:,mcycle)
      vr = vr(:,:,:,mcycle)
      call flux(vl, vr, f1d, ndir)
      F(:,:,:,mcycle,ndir) = f1d(:,:,:,:)

   end subroutine get_flux_ndir
end module timestep
