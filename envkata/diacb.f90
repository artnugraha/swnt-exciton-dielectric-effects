subroutine diacb(muq,q,nmk1,nmk2,nmwf2,nmwf3,&
wfcoe1,wfcoe2,wfcoe3,wfcoe4,vq2,con)
use common,only:nmmesh
use common,only:wfqp
use common,only:nc,nt
use common,only:wfmesh

use exkataif,only : selfcb

implicit none
!input
integer muq,q
complex*16 vq2(2,2)
type(nmmesh)::nmk1,nmk2
type(wfmesh)::nmwf2,nmwf3
type(wfqp)::wfcoe1,wfcoe2,wfcoe3,wfcoe4
!output
real*8 con
!local
!ft
complex*16,dimension(2,2)::vq
!self
real*8 selfenergy
!other
integer sn,sc
complex*16 vqt

vqt=0.0d0
do sn=1,2
do sc=1,2
vqt=vqt+vq2(sn,sc)
enddo
enddo

vqt=vqt/4.0d0


call selfcb(muq,q,nmk1,nmk2,nmwf2,nmwf3,&
wfcoe1,wfcoe2,wfcoe3,wfcoe4,selfenergy)


con=real(vqt)*selfenergy/(nc*nt)


end subroutine diacb



subroutine selfcb(muq,q,nmk1,nmk2,nmwf2,nmwf3,&
wfcoe1,wfcoe2,wfcoe3,wfcoe4,selfenergy)
use common,only:nc,nt
use common,only:nmmesh
use common,only:wfmesh
use common,only:wfqp
use common,only:nab,rab,hab,oab

use exkataif,only : cutoff
implicit none
!input
integer muq,q
type(nmmesh)::nmk1,nmk2
type(wfmesh)::nmwf2,nmwf3
type(wfqp)::wfcoe1,wfcoe2,wfcoe3,wfcoe4
!output
real*8 selfenergy
!local
!wavevector
integer muk,k
integer kc,kt
real*8,dimension(2)::kvec1,kvec2
!ec and ev
real*8 ec1,ev1,ec2,ev2
!wave function coefficient 
complex*16,dimension(2)::cc1,cv1,cc2,cv2
!tbtube
integer*4 nb,nd,nu,ns
parameter (nb=4,nd=3,nu=64,ns=2)
integer*4 nm
parameter (nm=nb*ns)
integer*4 evflag
parameter (evflag=1)
real*8 ener(nm)
complex*16 ampl(nm,nm)
!pol
real*8 f1d,f2d
complex*16 f1u,f2u
real*8 gv1,gv2,gc1,gc2,g1,g2,pol
!cutoff
!real*8,external::cutoff
!for check
integer mubt2,muup2,mubt3,muup3,kbt2,kup2,kbt3,kup3
logical cond
!other
integer i,ii


mubt2=nmwf2%mubt
muup2=nmwf2%muup
kbt2=nmwf2%kbt
kup2=nmwf2%kup

mubt3=nmwf3%mubt
muup3=nmwf3%muup
kbt3=nmwf3%kbt
kup3=nmwf3%kup

pol=0.0d0

do i=1,nmk1%nline
muk=nmk1%mu(i)
do k=nmk1%bt(i),nmk1%up(i),1

ev1=wfcoe1%ev(muk,k)
ec1=wfcoe1%ec(muk,k)
cv1(:)=wfcoe1%cv(muk,k,:)
cc1(:)=wfcoe1%cc(muk,k,:)


kc=muk+muq
kt=k+q

cond=kc.ge.mubt2.and.kc.le.muup2.and.kt.ge.kbt2.and.kt.le.kup2
if (.not.cond) then
write (*,100)
stop
endif

ev2=wfcoe3%ev(kc,kt)
ec2=wfcoe3%ec(kc,kt)
cv2(:)=wfcoe3%cv(kc,kt,:)
cc2(:)=wfcoe3%cc(kc,kt,:)


f1u=SUM(conjg(cv1)*cc2)

f1d=ec2-ev1

f2u=SUM(conjg(cc1)*cv2)
f2d=ec1-ev2

gv1=cutoff(ev1)
gv2=cutoff(ev2)
gc1=cutoff(ec1)
gc2=cutoff(ec2)
g1=gc2*gv1
g2=gc1*gv2

pol=pol+conjg(f1u)*f1u/f1d*g1+conjg(f2u)*f2u/f2d*g2

enddo
enddo


do i=1,nmk2%nline
muk=nmk2%mu(i)
do k=nmk2%bt(i),nmk2%up(i),1

ev1=wfcoe2%ev(muk,k)
ec1=wfcoe2%ec(muk,k)
cv1(:)=wfcoe2%cv(muk,k,:)
cc1(:)=wfcoe2%cc(muk,k,:)

kc=muk+muq
kt=k+q

cond=kc.ge.mubt3.and.kc.le.muup3.and.kt.ge.kbt3.and.kt.le.kup3
if (.not.cond) then
write (*,200)
stop
endif

ev2=wfcoe4%ev(kc,kt)
ec2=wfcoe4%ec(kc,kt)
cv2(:)=wfcoe4%cv(kc,kt,:)
cc2(:)=wfcoe4%cc(kc,kt,:)


f1u=SUM(conjg(cv1)*cc2)
f1d=ec2-ev1


f2u=SUM(conjg(cc1)*cv2)
f2d=ec1-ev2

gv1=cutoff(ev1)
gv2=cutoff(ev2)
gc1=cutoff(ec1)
gc2=cutoff(ec2)
g1=gc2*gv1
g2=gc1*gv2

pol=pol+conjg(f1u)*f1u/f1d*g1+conjg(f2u)*f2u/f2d*g2

enddo
enddo

selfenergy=2.0d0*pol


100 format (">>diacb>>wavevector range is not correct for K1")
200 format (">>diacb>>wavevector range is not correct for K2")

end subroutine selfcb
