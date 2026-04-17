subroutine kapfunc(ck,cex,aa,bb,cc,inpfile,outpfile)
! write energy and kappa with previously obtained function
! kapp = 1 + ck (p^aa/dt^bb . (1/lk)^cc - cex)

  implicit none

  ! general variables
  integer :: j, k, dummy
  integer, allocatable :: famNum(:) ! family
  
  ! variables for reading input of optimized kappa
  integer :: nlinei
  integer, allocatable :: n(:), m(:), p(:), i(:)
  real(8), allocatable :: dt(:), theta(:), Eexp(:), lk(:), kapp(:)
  character(30), intent(in) :: inpfile
  character(30) :: outpfile

  ! output
  real(8), allocatable :: Ecal(:), diffE(:)

  ! variables for data processing
  real(8), intent(in) :: ck, cex, aa, bb, cc 

  ! variables for reading sx-eXX.dat:
  integer :: nliner
  integer, allocatable :: nr(:), mr(:)
  real(8), allocatable :: dr(:), thetar(:), Er(:), lkr(:), kappr(:)
  character(25) :: filetoread
  character(15), parameter :: root = '~/repo/'

  ! inpfile format: n m dt theta E lk optkapp p i
  open(1,file=trim(inpfile),status='old')
  nlinei = 0
  do
     read(1,*,end=10), dummy
     nlinei = nlinei+1
  end do
10 rewind(1)
  allocate(n(nlinei))
  allocate(m(nlinei))
  allocate(dt(nlinei))
  allocate(theta(nlinei))
  allocate(Eexp(nlinei))
  allocate(Ecal(nlinei))
  allocate(lk(nlinei))
  allocate(kapp(nlinei))
  allocate(p(nlinei))
  allocate(i(nlinei))
  allocate(diffE(nlinei))
  allocate(famNum(nlinei))

  open(3,file=outpfile)
  do j = 1, nlinei
     read(1,*), n(j), m(j), dt(j), theta(j), Eexp(j),&
          & lk(j), dummy, p(j), i(j)  
     kapp(j) = 1 + ck*((p(j)*aa)/(dt(j)**bb) * (1./(lk(j)**cc)) - cex)
     open(4,file='tmp.dat')
     write(4,'(f6.1)'), kapp(j)
     rewind(4); read(4,'(f6.1)'), kapp(j)
     close(4,status='delete')
     famNum(j) = 2*n(j) + m(j)
     if (mod(famNum(j),3) == 0) then
        filetoread = trim(root)//'s0-e'//char(48+i(j))//char(48+i(j))//'.dat'
     else if (mod(famNum(j),3) == 1) then
        filetoread = trim(root)//'s1-e'//char(48+i(j))//char(48+i(j))//'.dat'
     else if (mod(famNum(j),3) == 2) then
        filetoread = trim(root)//'s2-e'//char(48+i(j))//char(48+i(j))//'.dat'
     else
        print*, 'incorrect family number!'
     end if
     
     ! file to read: n, m, dt, theta, E(cal), lk, kappa
     open(2,file=filetoread,status='old')
     nliner = 0
     do
        read(2,*,end=20), dummy
        nliner = nliner+1
     end do
     print*, nliner
20   rewind(2)
     allocate(nr(nliner))
     allocate(mr(nliner))
     allocate(dr(nliner))
     allocate(thetar(nliner))
     allocate(Er(nliner))
     allocate(lkr(nliner))
     allocate(kappr(nliner))

     do k = 1, nliner
        read(2,*), nr(k), mr(k), dr(k), thetar(k), Er(k), &
             & lkr(k), kappr(k)
        if ((nr(k) == n(j) .and. mr(k) == m(j)) .and.&
             &kappr(k) == kapp(j)) then
           Ecal(j) = Er(k)
           
        end if
     end do
     diffE(j) = Eexp(j) - Ecal(j)
     if (abs(int(1000*diffE(j))) <= 70) then
        write(3,'(2i6,5f10.3,f7.1,3i6)'), n(j), m(j), dt(j), theta(j),&
             & Ecal(j), Eexp(j), lk(j), kapp(j), int(1000.*diffE(j)),&
             & p(j), i(j)
     end if
     deallocate(nr,mr,dr,thetar,Er,lkr,kappr)
     close(2)
  end do

  close(1)
  close(3)

  deallocate(n,m,dt,theta,Eexp,Ecal,lk,kapp,p,i,diffE)

end subroutine
