program linkapp
! find "basic" constant ckapp and cex in the following equation:
! kapp - 1 = ck (p^aa/dt^bb . (1/lk)^cc - cex)
!
! equivalently,
! define y = kapp-1 and x = p^aa/dt^bb.(1/lk)^cc, A = -ck*cex, B = ck
!
! then we simply express yi = A + B*xi
!
! after that the ck and cex are varied to find minimum R based on
! experimental data
!
!======================================================================!
! VARIABLES DECLARATION

  implicit none

  ! general variables
  integer :: j, k, l, dummy
  integer, parameter :: mesh = 100
  
  ! variables for reading input of optimized kappa
  integer :: nlinei
  integer, allocatable :: n(:), m(:), p(:), i(:)
  real(8), allocatable :: dt(:), Eexp(:), Ecal(:), lk(:), optkapp(:)
  character(30) :: inputfile

  ! variables for data processing
  real(8), allocatable :: x(:), y(:)
  real(8) :: sumx, sumy, sumx2, sumy2, sumxy, A, B
  real(8) :: ssxx, ssyy, ssxy, sstd
  real(8), parameter :: aa = 0.8, bb = 1.6, cc = 0.4

  ! variables for intermediate output
  real(8) :: ck, cex, Rsqr, Bstderr

  ! variables for intermediate processing
  real(8) :: ckmin, ckmax, cexmin, cexmax, Reval

!======================================================================!
! INPUT: e.g. opt-hata.dat

  !inputfile = "opt-hata.dat"
  inputfile = "hata-opt.dat"

  ! inputfile format: n m dt Eexp lk optkapp p
  open(1,file=trim(inputfile),status="old")
  nlinei = 0
  do
     read(1,*,end=10), dummy
     nlinei = nlinei+1
  end do
10 rewind(1)
  allocate(n(nlinei))
  allocate(m(nlinei))
  allocate(dt(nlinei))
  allocate(Eexp(nlinei))
  allocate(Ecal(nlinei))
  allocate(lk(nlinei))
  allocate(optkapp(nlinei))
  allocate(p(nlinei))
  allocate(i(nlinei))
  allocate(x(nlinei))
  allocate(y(nlinei))
  
!======================================================================!
! DATA PROCESSING

  ! finding basic ck and cex

  ! initialization
  sumx = 0.d0; sumy = 0.d0; sumx2 = 0.d0; sumxy = 0.d0
  do j = 1, nlinei
     read(1,*), n(j), m(j), dt(j), Eexp(j), lk(j),&
          & optkapp(j), p(j), i(j)
     x(j) = (p(j)*aa)/(dt(j)**bb) * (1./(lk(j)**cc))
     y(j) = optkapp(j) - 1.d0
     sumx = sumx + x(j)
     sumy = sumy + y(j)
     sumx2 = sumx2 + x(j)**2
     sumy2 = sumy2 + y(j)**2
     sumxy = sumxy + x(j)*y(j)
  end do
  close(1)

  ! linear regression formula
  A = (sumx2*sumy - sumx*sumxy) / (nlinei*sumx2 - (sumx)**2)
  B = (nlinei*sumxy - sumx*sumy) / (nlinei*sumx2 - (sumx)**2)

  ssxx = sumx2 - ( (sumx**2) / nlinei )
  ssyy = sumy2 - ( (sumy**2) / nlinei )
  ssxy = sumxy - ( (sumx*sumy/ nlinei))
  sstd = sqrt((ssyy-(ssxy**2/ssxx))/dble(nlinei-2))

  Bstderr = sstd / sqrt(ssxx)

  Rsqr = (sumxy - sumx*sumy/nlinei)**2 / &
       & ((sumx2 - sumx**2/nlinei)*(sumy2 - sumy**2/nlinei))

!----------------------------------------------------------------------!

  ! write ck and cex
  ck = B
  cex = - A / B

  print"(a,f10.3)", "Rsqr = ", Rsqr
  print"(a,2f7.3)", "# Coefficient ck and cex = ", ck, cex

!----------------------------------------------------------------------!

  ! define ckmin -> ckmax, cexmin -> cexmax
  cexmin = 1.d0
  cexmax = 2.d0
  ckmin = ck - 0.1d0
  ckmax = ck + 0.1d0
  
  ! now ck and cex is no longer previous constant, but varied

  ! initialization
  cex = cexmin
  ck = ckmin
  open(3,file="meshhataorig.dat")
  do j = 1, mesh+1
     do k = 1, mesh+1
        Reval = 0.d0
        do l = 1, nlinei
           y(l) = ck*(x(l) - cex)
           optkapp(l) = y(l) + 1.d0
           open(2,file="tmporig.dat")
           write(2,"(f7.1)"), optkapp(l)
           rewind(2)
           read(2,"(f7.1)"), optkapp(l)
           close(2,status="delete")
           call findeii(n(l),m(l),p(l),i(l),optkapp(l),Ecal(l))
           print"(f7.1,f7.3)", optkapp(l),Ecal(l)
           Reval = Reval + (1000*(Eexp(l) - Ecal(l))/10)**2 ! 10 = errorbar
        end do
        write(3,"(2f7.3,f12.3)"), cex, ck, Reval
        ck = ckmin + k*(ckmax-ckmin)/mesh
     end do
     cex = cexmin + j*(cexmax-cexmin)/mesh
  end do
  close(3)

  deallocate(n,m,dt,Eexp,Ecal,lk,optkapp,p,i,x,y)

!======================================================================!

end program
