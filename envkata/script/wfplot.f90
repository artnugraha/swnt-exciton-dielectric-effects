
subroutine wfplot(otype,tubetype,cline,vline,state,nd,vec,kcv,fcv)

  use common,only:ns,nu
  use common,only:nab
  use common,only:bound
  use common,only:rabc,rabn
  use common,only:xyzc,xyzn
  use common,only:dt

  use common,only : nn, mm

  implicit none

  ! input
  character(len=2) :: otype, tubetype
  integer :: cline, vline, state
  complex(8) :: vec(bound)
  integer :: nd
  real(8) :: kcv(bound,2,2,2)
  complex(8) :: fcv(bound,2,2,2)

  ! local
  integer :: nsite,sc,sn,un,u
  complex(8) :: iu
  complex(8) :: zkc1,zkv1,zkc2,zkv2
  complex(8) :: ph1,ph2,base1,base2,base
  complex(8) :: wf
  integer :: i

  ! constants
  real(8),parameter :: pi=3.1415926d0,sq3=1.732d0

  real(8),dimension(2)::kc1,kv1,kc2,kv2,lth,lte,xyh,xye,xyei,lteimg,lthimg
  integer num(-100:100),nmax,dn,up,ntest
  real(8) :: xref(-100:100),xi
  integer,parameter::natom=5000
  real(8) :: xyef(-100:100,natom,2),ltef(-100:100,natom,2),xyzf(-100:100,natom,3)
  real(8),parameter :: eps=0.1d0
  character(len=12) :: file
  integer unum,n1,n2
  real(8) :: xsite,yadd
  real(8) :: acc,auc

  !integer u,i1,nadd
  integer i1,nadd

  parameter (acc=0.142d0,auc=sq3*acc)
  integer unum1,ii,jj
  real(8) :: th,wfr,norm,x,y

  integer,parameter :: AsciiZero = 48
  character(len=100) :: fn100, fn101, fn102

  
  fn100 = './exwf/axiswf-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
       &//char(state-1+AsciiZero)//'.dat'
  fn101 = './exwf/cirwf-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
       &//char(state-1+AsciiZero)//'.dat'
  fn102 = './exwf/cirwffit-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
       &//char(state-1+AsciiZero)//'.dat'

  open(100,file=fn100)
  open(200,file=fn101)
  open(300,file=fn102)

  sc=1
  sn=1

  nadd=nab(sc,sc)+1
  nab(sc,sc)=nab(sc,sc)+1
  rabn(:,nadd,sc,sc)=rabc(:,sc)
  xyzn(:,nadd,sc,sc)=xyzc(:,sc)

  iu=cmplx(0.0d0,1.0d0)

  nsite = 0
  do un = 1, nab(sn,sc)
     nsite = nsite + 1
  enddo

  !dn=-7
  !up=8

  dn=-31
  up=32

  xi=-0.25
  do u=dn,up
     xref(u)=xi+0.5d0*u
  enddo

  num = 0
  un1doloop:do un = 1, nab(sn,sc)
     xyei(:) = (/ rabn(1,un,sn,sc)*dt*0.5d0,rabn(2,un,sn,sc) /)
     do u = dn, up
        if( abs(xyei(1)-xref(u)).lt.eps ) then
           num(u) = num(u) + 1
           if (num(u).gt.natom) then
              write (*,500) 
           endif
           xyef(u,num(u),:) = xyei(:)
           ltef(u,num(u),:) = rabn(:,un,sn,sc)
           xyzf(u,num(u),:) = xyzn(:,un,sn,sc)
        endif
     enddo
  enddo un1doloop

  ntest = 0
  do u = dn, up
     ntest = ntest + num(u)
  enddo

  if (ntest.ne.nsite) then
     write (*,*) ntest, nsite
     write (*,400) 
     stop
  endif

  if (sn.eq.1) then
     xyh(:)=(/rabc(1,sc)*dt*0.5d0,rabc(2,sc)/)
     lth(:)=rabc(:,sc)
  elseif (sn.eq.2) then
     xyh(:)=(/rabn(1,1576,2,sc)*dt*0.5d0,rabn(2,1576,2,sc)/)
     lth(:)=rabn(:,1576,2,sc)
  endif


  !for axis direction
  unum = 0
  uloop1:do u = dn, up
     unum = unum + 1
     xsite = unum*2.0d0
     yadd = -0.866*(mod(unum,2)-1)

     unum1 = 0
     un2loop1:do un = 1, num(u)
        xye(:) = xyef(u,un,:)
        lte(:) = ltef(u,un,:)

        if (abs(xye(2)).gt.100.0d0) cycle
        if (abs(lte(1)-lth(1)).gt.pi/16.0d0) cycle

        wf = dcmplx(0.0d0,0.0d0)

        do i = 1, nd
           kv1(:) = kcv(i,1,1,:)
           kc1(:) = kcv(i,2,1,:)
           zkv1 = fcv(i,1,1,sn)
           zkc1 = fcv(i,2,1,sn)
           ph1 = exp(iu*(SUM(kc1*lte)-SUM(kv1*lth)))
           kv2(:) = kcv(i,1,2,:)
           kc2(:) = kcv(i,2,2,:)
           zkv2 = fcv(i,1,2,sn)
           zkc2 = fcv(i,2,2,sn)
           ph2 = exp(iu*(SUM(kc2*lte)-SUM(kv2*lth)))
           base1 = zkc1*conjg(zkv1)
           base2 = zkc2*conjg(zkv2)
           select case(otype)
           case('ge')
              base = base1*ph1 - base2*ph2
           case('go')
              base = base1*ph1 + base2*ph2
           case ('k1')
              base = base1*ph1
           case('k2')
              base = base1*ph1
           end select

           wf = wf + vec(i)*base
        enddo

        write (100,'(10f12.4)') (lte(2)+yadd)*auc,real(wf)

     enddo un2loop1

     if (unum1.ne.0) then
        write (100,*) '&'
     endif


  enddo uloop1


  !for circumference direction

  unum=0
  uloop2:do u=dn,up
     unum=unum+1
     xsite=unum*2.0d0
     yadd=-0.866*(mod(unum,2)-1)
     un2loop2:do un=1,num(u)
        xye(:)=xyef(u,un,:)
        lte(:)=ltef(u,un,:)

        if (abs(xye(2)-xyh(2)).gt.1.0d0/(3.0d0*sqrt(3.0d0))) cycle

        wf=dcmplx(0.0d0,0.0d0)
        do i=1,nd
           kv1(:)=kcv(i,1,1,:)
           kc1(:)=kcv(i,2,1,:)
           zkv1=fcv(i,1,1,sn)
           zkc1=fcv(i,2,1,sn)
           ph1=exp(iu*(SUM(kc1*lte)-SUM(kv1*lth)))
           kv2(:)=kcv(i,1,2,:)
           kc2(:)=kcv(i,2,2,:)
           zkv2=fcv(i,1,2,sn)
           zkc2=fcv(i,2,2,sn)
           ph2=exp(iu*(SUM(kc2*lte)-SUM(kv2*lth)))
           base1=zkc1*conjg(zkv1)
           base2=zkc2*conjg(zkv2)
           select case(otype)
           case('ge')
              base=base1*ph1-base2*ph2
           case('go')
              base=base1*ph1+base2*ph2
           case ('k1')
              base=base1*ph1
           case('k2')
              base=base1*ph1
           end select
           wf=wf+vec(i)*base

        enddo


        write (200,'(10f12.4)') lte(1),real(wf),imag(wf)


     enddo un2loop2
  enddo uloop2




  !fitting for circumference direction

  lteimg(2)=lth(2)

  do ii=0,200

     lteimg(1)=2.0d0*pi/200.0d0*ii

     wf=dcmplx(0.0d0,0.0d0)
     do i=1,nd
        kv1(:)=kcv(i,1,1,:)
        kc1(:)=kcv(i,2,1,:)
        zkv1=fcv(i,1,1,sn)
        zkc1=fcv(i,2,1,sn)
        ph1=exp(iu*(SUM(kc1*lteimg)-SUM(kv1*lth)))
        kv2(:)=kcv(i,1,2,:)
        kc2(:)=kcv(i,2,2,:)
        zkv2=fcv(i,1,2,sn)
        zkc2=fcv(i,2,2,sn)
        ph2=exp(iu*(SUM(kc2*lteimg)-SUM(kv2*lth)))
        base1=zkc1*conjg(zkv1)
        base2=zkc2*conjg(zkv2)
        select case(otype)
        case('ge')
           base=base1*ph1-base2*ph2
        case('go')
           base=base1*ph1+base2*ph2
        case ('k1')
           base=base1*ph1
        case('k2')
           base=base1*ph1
        end select
        wf=wf+vec(i)*base
     enddo

     write (300,'(10f12.4)') lteimg(1),real(wf),imag(wf)

  enddo


  !dipole for circumference direction

  lteimg(2)=lth(2)
  lthimg(2)=lth(2)

  th=0.0d0
  norm=0.0d0

  do ii=0,200

     lthimg(1)=2.0d0*pi/200.0d0*ii
     !lthimg(1)=0.0d0

     do jj=0,200

        lteimg(1)=2.0d0*pi/200.0d0*jj

        wf=dcmplx(0.0d0,0.0d0)
        do i=1,nd
           kv1(:)=kcv(i,1,1,:)
           kc1(:)=kcv(i,2,1,:)
           zkv1=fcv(i,1,1,sn)
           zkc1=fcv(i,2,1,sn)
           ph1=exp(iu*(SUM(kc1*lteimg)-SUM(kv1*lthimg)))
           kv2(:)=kcv(i,1,2,:)
           kc2(:)=kcv(i,2,2,:)
           zkv2=fcv(i,1,2,sn)
           zkc2=fcv(i,2,2,sn)
           ph2=exp(iu*(SUM(kc2*lteimg)-SUM(kv2*lthimg)))
           base1=zkc1*conjg(zkv1)
           base2=zkc2*conjg(zkv2)
           select case(otype)
           case('ge')
              base=base1*ph1-base2*ph2
           case('go')
              base=base1*ph1+base2*ph2
           case ('k1')
              base=base1*ph1
           case('k2')
              base=base1*ph1
           end select
           wf=wf+vec(i)*base
        enddo

        wfr=real(wf)

        norm=norm+wfr**2

        th=th+(lteimg(1)-lthimg(1))*wfr**2
        x=x+dt*0.5d0*cos(lteimg(1))*wfr**2
        y=y+dt*0.5d0*sin(lteimg(1))*wfr**2


     enddo

  enddo

  write (*,*) 'ok',x/norm,y/norm

  close(100)
  close(200)
  close(300)

500 format(" >> wfplot >> natom too small >> ")
400 format(" >> wfplot >> number of site error >> ")

  return
end subroutine wfplot
