!program testeii
!
!  implicit none
!
!  integer :: nn, mm, pp, ii
!  real(8) :: kappa, Ethe
!
!  nn = 10
!  mm = 7
!  ii = 2
!  pp = 3
!  kappa = 3.5
!
!  call findeii(nn,mm,pp,ii,kappa,Ethe)
!  print"(a,f7.3)", "Eth = ", Ethe
!
!end program testeii

subroutine findeii(n,m,p,i,kapp,Eth)
! find Eii from kappa

!======================================================================!
! VARIABLES DECLARATION

  implicit none

  ! general variables
  integer :: k, dummy
  integer, intent(in) :: n, m, p, i
  real(8), intent(in) :: kapp

  ! variables for reading sx-eXX.dat:
  integer :: nliner, famNum
  integer, allocatable :: nr(:), mr(:)
  real(8), allocatable :: dr(:), thetar(:), Er(:), lkr(:), kappr(:)
  character(25) :: filetoread
  character(15), parameter :: root = "genkap/"

  ! output
  real(8), intent(out) :: Eth

!======================================================================!

! READING DATABASE

  famNum = 2*n + m
  if (mod(famNum,3) == 0) then
     filetoread = trim(root)//'s0-e'//char(48+i)//char(48+i)//'.dat'
  else if (mod(famNum,3) == 1) then
     filetoread = trim(root)//'s1-e'//char(48+i)//char(48+i)//'.dat'
  else if (mod(famNum,3) == 2) then
     filetoread = trim(root)//'s2-e'//char(48+i)//char(48+i)//'.dat'
  else
     print*, 'incorrect family number!'
  end if

  ! file to read: n, m, dt, theta, E(cal), lk, kappa
  open(2,file=trim(filetoread),status="old")
  nliner = 0
  do
     read(2,*,end=20), dummy
     nliner = nliner+1
  end do
  print*, nliner
20 rewind(2)
  allocate(nr(nliner))
  allocate(mr(nliner))
  allocate(dr(nliner))
  allocate(thetar(nliner))
  allocate(Er(nliner))
  allocate(lkr(nliner))
  allocate(kappr(nliner))

!----------------------------------------------------------------------!
! set the condition to find Eii from database

  do k = 1, nliner
     read(2,*), nr(k), mr(k), dr(k), thetar(k), Er(k), &
          & lkr(k), kappr(k)
     if ((n == nr(k) .and. m == mr(k))&
          & .and. nint(10.*kapp) == nint(10.*kappr(k))) then
        Eth = Er(k)
     end if
  end do
!----------------------------------------------------------------------!

  deallocate(nr,mr,dr,thetar,Er,lkr,kappr)
  close(2)

!======================================================================!

end subroutine findeii
