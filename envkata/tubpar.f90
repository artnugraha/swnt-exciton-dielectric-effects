!>
!! \file tubpar.f90
!! \brief cntgrlib (Carbon Nanotube and Graphite Library)
!<

!>
!! \ingroup tbdftsnt
!! \brief structural parameters of nanotubes rolled up from graphene
!!
!! \param[in] n nanotube structural indices
!! \param[in] m nanotube structural indices
!! \param[in] l nanotube length [nm]
!!
!! \param[out] t  nanotube unit cell length [nm]
!! \param[out] d  nanotube diameter [nm]
!! \param[out] theta  nanotube chiral angle [radian]
!! \param[out] tab  translational distance between A & B atoms [t]
!! \param[out] phiab  circumferential distance between A & B atoms [radian]
!! \param[out] dab  radial distance between A & B atoms [d]
!! \param[out] a1t  translational component of graphene unit vector a1 [t]
!! \param[out] a1phi  circumferential component of a1 vector [radian]
!! \param[out] a2t  translational component of graphene unit vector a2 [t]
!! \param[out] a2phi  circumferential component of a2 vector [radian]
!! \param[out] nc  number of cutting lines in circumferential direction
!! \param[out] nt  number of mesh points in translational direction
!! \param[out] kcx,kcy  circumferential quantization wavevector [1/auc]
!! \param[out] ktx,kty  translational quantization wavevector [1/auc]
!!
!! \notice
!! - tab=(tb-ta)/t & phiab=phib-phia & dab=(db-da)/d & d=(db+da)/2
!!   parameters (t,d,tab,phiab) need to be optimized.
!! - parameters (dab) may also be optimized
!! - parameters (a1t,a1phi,a2t,a2phi) are fixed
!!
!! \author Georgii Samsonidze <gsm@mgm.mit.edu>
!!
!! \date last modified 2004-07-02
!! \date created 2004-05-11
!<
subroutine tubpar(n,m,l,t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nc,nt,kcx,kcy,ktx,kty)

  implicit none

  integer*4 n,m
  real*8 l
  real*8 t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi
  integer*4 nc,nt
  real*8 kcx,kcy,ktx,kty
  real*8 pi,sr3
  parameter (pi=+0.3141592653589793d+01,sr3=+0.1732050807568877d+01)
  real*8 acc,auc
  parameter (acc=+0.142d+00,auc=acc*sr3)
  real*8 eps
  parameter (eps=+0.1d-08)
  integer*4 dr,t1,t2
  real*8 r

  if ((n.eq.0.and.m.eq.0).or.l.le.eps) then
     write(*,100)
     stop
  endif
  call gcd (2*n+m,n+2*m,dr)
  t=auc*sr3*dsqrt(dfloat(n*n+n*m+m*m))/dr
  d=auc*dsqrt(dfloat(n*n+n*m+m*m))/pi
  theta=datan(sr3*dfloat(m)/dfloat(2*n+m))
  r=d/2.0d0
  tab=dsin(theta-pi/6.0d0)*acc/t
  phiab=dcos(theta-pi/6.0d0)*acc/r
  dab=0.0d0
  a1t=dsin(theta)*auc/t
  a1phi=dcos(theta)*auc/r
  a2t=dsin(theta-pi/3.0d0)*auc/t
  a2phi=dcos(theta-pi/3.0d0)*auc/r
  nc=2*(n*n+n*m+m*m)/dr
  nt=2*(idint(l/t)/2)
  t1=dfloat((2*m+n)/dr)
  t2=-dfloat((2*n+m)/dr)
  kcx=2.0d0*pi/sr3*(-t2+t1)/dfloat(nc)
  kcy=2.0d0*pi*(-t2-t1)/dfloat(nc)
  ktx=2.0d0*pi/sr3*dfloat(m-n)/dfloat(nc)
  kty=2.0d0*pi*dfloat(m+n)/dfloat(nc)
  ktx=ktx/dfloat(nt)
  kty=kty/dfloat(nt)

100 format(" >> cntgrlib >> tubpar >> out of range >> ")
end subroutine tubpar
