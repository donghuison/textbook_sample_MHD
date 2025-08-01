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
      call step_split
   end subroutine step_full


#define STEP_SPLIT_1D(n) step_split_rk2_1d(n)

   subroutine step_split

      use grid, only: V, W, F, Dtime
      use flux_eos, only : source_b, v2u, u2v

      call step_split_2d
      call v2u(V, W)
      call source_b(F, W, Dtime)
      call u2v(W, V)

   end subroutine step_split

   subroutine step_split_2d
      integer :: n
      do n = MX, MY
         call STEP_SPLIT_1D( n )
      end do
   end subroutine step_split_2d

   subroutine step_split_rk2_1d(ndir)
      use grid
      use flux_eos
      use boundary

      integer,intent(IN) :: ndir

      call v2u(V, U)
      call boundary_fix(V)
      call get_flux_ndir(ndir)
      W = U
      call w_update_ndir(Dtime, ndir)
      call u2v(W, V)
      call boundary_fix(V)
      call get_flux_ndir(ndir)
      W = (U + W)*0.5d0
      call w_update_ndir(Dtime*0.5d0, ndir)
      call u2v(W, V)
      call boundary_fix(V)
   end subroutine step_split_rk2_1d

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

      if (b_muscl) then
         do m = MMIN, MMAX
            do k = KMIN-ko, KMAX
               do j = JMIN-jo, JMAX
                  do i = IMIN-io, IMAX
                     vl(i,j,k,m) = vl(i,j,k,m) &
                        + (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i,j,k,m)-V(i-io,j-jo,k-ko,m)))*0.5d0
                     vr(i,j,k,m) = vr(i,j,k,m) &
                        - (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i+i2,j+j2,k+k2,m)-V(i+io,j+jo,k+ko,m)))*0.5d0
                  enddo
               enddo
            enddo
         enddo
      endif
      mcycle = cyclecomp( ndir )
      vl = vl(:,:,:,mcycle)
      vr = vr(:,:,:,mcycle)
      call flux(vl, vr, f1d, ndir)
      F(:,:,:,mcycle,ndir) = f1d(:,:,:,:)

   end subroutine get_flux_ndir
end module timestep
