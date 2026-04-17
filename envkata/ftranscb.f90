subroutine ftranscb(mu,q,type,vq)
use common,only:nc,nt
use common,only:u0,v0
use common,only:ns,nabf,rabc,rabn,xyzc,xyzn
implicit none
!  input
integer*4 type
real*8 mu,q
!  output
complex*16 vq(ns,ns)
! local
integer*4 sc,sn,un
real*8 xc,yc,zc,xn,yn,zn,x,y,z,r,vr
!wavevector
real*8 fkt,fkphi
!parameter
real*8 pi,sr3
parameter (pi=+0.3141592653589793d+01,sr3=+0.1732050807568877d+01)
real*8 acc,auc
parameter (acc=+0.142d+00,auc=acc*sr3)
!phase factor
real*8 rab(2)
complex*16 pf



fkphi=mu
fkt=q/dfloat(nt)*2.0d0*pi

vq=dcmplx(0.0d0,0.0d0)

do sc=1,ns
xc=xyzc(1,sc)
yc=xyzc(2,sc)
zc=xyzc(3,sc)
do sn=1,ns

do un=1,nabf(sn,sc)
xn=xyzn(1,un,sn,sc)
yn=xyzn(2,un,sn,sc)
zn=xyzn(3,un,sn,sc)
x=xn-xc
y=yn-yc
z=zn-zc
r=sqrt(x*x+y*y+z*z)
!
!r=r/acc
!
r=r/auc
!
vr=v0/r
rab(:)=rabn(:,un,sn,sc)-rabc(:,sc)
pf=cdexp(dcmplx(0.0d0,fkphi*rab(1)+fkt*rab(2)))
vq(sn,sc)=vq(sn,sc)+pf*vr
enddo

if (sn.eq.sc.and.type.eq.2) then
vq(sn,sc)=vq(sn,sc)+u0
endif


enddo
enddo

return
end subroutine ftranscb
