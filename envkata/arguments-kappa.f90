
subroutine GetArguments(type,subband,kappa)

  implicit none

  character(len=2),intent(out) :: type
  ! nanotube type
  ! s0: metal
  ! s1: Type I
  ! s2: Type II

  integer,intent(out) :: subband
  ! subband index: 1, 2, 3...

  real(8),intent(out) :: kappa
  ! dielectric constant

  ! arguments
  integer,parameter :: BUFSIZE = 100
  character(len=BUFSIZE) :: buf
  integer,parameter :: an1 = 1, an2 = 2, an3 = 3
  integer :: in_kappa

  call getarg(an1,buf)
  read(buf,FMT='(A2)') type

  call getarg(an2,buf)
  read(buf,FMT='(I3)') subband

  call getarg(an3,buf)
  read(buf,FMT='(I4)') in_kappa
  kappa = dble(in_kappa)*0.1D0

end subroutine GetArguments



subroutine SRArguments(my_rank,p,type,subband,kappa)

  implicit none

  ! header of MPI
  include "mpif.h"

  ! input
  integer,intent(in) :: my_rank, p

  ! output
  character(len=2),intent(out) :: type
  integer,intent(out) :: subband
  real(8),intent(out) :: kappa

  ! for MPI
  integer :: status(MPI_STATUS_SIZE), ierr
  ! local
  integer,parameter :: cpuzero = 0, count1 = 2, count2 = 1, count3 = 1, tag1 = 1, tag2 = 2, tag3 = 3
  integer :: i

  if( my_rank /= cpuzero ) then
     ! recieve arguments from CPU zero
     call MPI_Recv( type, count1, MPI_CHARACTER, cpuzero, tag1, MPI_COMM_WORLD, status, ierr )
     call MPI_Recv( subband, count2, MPI_INTEGER, cpuzero, tag2, MPI_COMM_WORLD, status, ierr )
     call MPI_Recv( kappa, count3, MPI_REAL8, cpuzero, tag3, MPI_COMM_WORLD, status, ierr )
  else if( my_rank == cpuzero ) then
     ! get arguments on CPU zero
     call GetArguments(type,subband,kappa)
     ! send arguments to other CPU
     do i = 1, p-1
        call MPI_Send( type, count1, MPI_CHARACTER, i, tag1, MPI_COMM_WORLD, ierr )
        call MPI_Send( subband, count2, MPI_INTEGER, i, tag2, MPI_COMM_WORLD, ierr )
        call MPI_Send( kappa, count3, MPI_REAL8, i, tag3, MPI_COMM_WORLD, ierr )
     end do
  end if

end subroutine SRArguments
