
include 'kappa.f90'

program main

  implicit none

  character(len=2) :: tubetype
  integer :: cline
  real(8) :: dt
  real(8) :: kappa
  integer :: p

  tubetype = 's1'
  cline = 2

  read(*,*) dt
  call SelectKappa(tubetype,cline,dt,kappa)

end program main
