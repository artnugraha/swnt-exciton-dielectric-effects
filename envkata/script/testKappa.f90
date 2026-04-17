
include 'kappa.f90'

program main

  implicit none

  real(8) :: dt, d
  real(8) :: pdt, kappa
  integer,parameter :: maxi = 100
  real(8),parameter :: mindt = 0.5D0, ddt = 0.1D0
  real(8),parameter :: KappaEnv= 3.0D0
  integer :: i

  do i = 0, maxi
     dt = mindt + i*ddt
     pdt = dble(1) / dt
     call KappaFunction(pdt,KappaEnv,Kappa)
     write(*,FMT='(4F12.5)') dt, pdt, KappaEnv, Kappa
  end do
  
end program main
