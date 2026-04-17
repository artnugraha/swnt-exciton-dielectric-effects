!>
!! \file num2str.f90
!! \brief change numbers to characters
!<

!>
!! \brief change numbers to characters
!!
!! \param[in] nn Chiral index n of (n,m)
!! \param[in] mm Chiral index m of (n,m)
!! \param[out] num2str character NNMM
!!
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
character(4) function num2str(nn,mm)

  integer,intent(in) :: nn, mm

  integer :: k
  character(len=1) :: n1, n2, m1, m2

  ! Ascii character number of Zero. See "man ascii".
  integer,parameter :: AZ = 48

  if( (nn<0).or.(99<nn) ) stop 'num2str expects n=0..99'
  if( (mm<0).or.(nn<mm) ) stop 'num2str expects m=0..n'

  k = nn / 10
  n1 = char(AZ+k)

  k = mod(nn,10)
  n2 = char(AZ+k)

  k = mm / 10
  m1 = char(AZ+k)

  k = mod(mm,10)
  m2 = char(AZ+k)

  num2str = n1//n2//m1//m2

end function num2str
