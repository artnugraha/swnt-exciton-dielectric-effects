!
!== Description
! 1. Calculate pi band energy and wave function coefficient based on STB
!
!== History
! This subroutine name EG3 has been used anywhere in RS Group
! for a histrical reason. Probably you find the similar name 
! subroutine such as EG or EG2...
!
!== References
! 1. Chapeter 2, R. Saito et al., Physical Properties of Carbon Nanotubes.
!
!== Contact List
! 1. Jie Jiang (mailto:jiang@chips.ncsu.edu)
! 2. Kentaro Sato (mailto:kentaro@flex.phys.tohoku.ac.jp)
!
function eg3(x,y)
  !== calculate energy and wave function coefficient of pi band
  !
  ! * K1 and K2 in units 1/a
  ! * mu 0..N-1 (cutting line index. different from previous def.)
  ! * k -1/2...1/2 (1D k along cutting line)

  use exkataparameter,only : I_U, SR3, TPPP, SPPP, EGresult

  implicit none

  real(8),intent(in) :: x, y
  ! wave vector of electron. unit is 1/a. 
  ! a is lattice constant. 

  type(EGresult) :: eg3
  ! pi band energy and wave function coefficient.
  ! See object composition for tight binding, EGresult.

  real(8) :: w
  complex(8) :: f, ca(2), cb(2)

  f = exp(i_u*x/SR3) + 2.0D0*exp(-i_u*x/(2*SR3))*cos(y/2.0D0)
  w = zabs(f)

  ca(1:2) = (/ sqrt(0.5d0/(1.0D0+sppp*w))*sqrt(f/w), sqrt(0.5d0/(1.0D0-sppp*w))*sqrt(f/w) /)
  cb(1:2) = (/ sqrt(0.5d0/(1.0D0+sppp*w))*sqrt(Conjg(f)/w), -sqrt(0.5d0/(1.0D0-sppp*w))*sqrt(Conjg(f)/w) /)
  ! TPPP is negative thus the plus sign appears for w1 which is valence band energy.

  eg3%w1 =  tppp*w / (1.0D0+sppp*w)
  eg3%w2 = -tppp*w / (1.0D0-sppp*w)
  eg3%ca(:) = ca(:)
  eg3%cb(:) = cb(:)

end function eg3
