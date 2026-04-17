!>
!! \file
!<

!>
!! \ingroup tbdftsnt
!! \brief generate nanotube coordinates
!<
subroutine tbcoord

  use common,only:ns,nu,nabf,rabc,rabn,xyzc,xyzn
  use common,only:tt,dt,tab,phiab,dab,a1t,a1phi,a2t,a2phi

  use exkataparameter,only : PI, SR3, ACC, AUC, cutlength

  implicit none

  ! local
  real(8) :: t,d
  real(8),parameter :: eps = 0.1D-8
  integer :: sc,sn,un,unmax,un1,un2
  real(8) :: cylcut,fsc,fsn,fun1,fun2
  real(8) :: tc,phic,dc,tn,phin,dn,xc,yc,zc,xn,yn,zn,x,y,z,r

  t = tt
  d = dt

  if( d.gt.cutlength ) then
     cylcut = d*dasin(cutlength/d)
  elseif( d.lt.cutlength ) then
     cylcut = dsqrt((d*pi/2.0d0)**2+(cutlength**2-d**2))
  else
     cylcut = d*pi/2.0d0
  end if

  unmax = idint((cylcut*2.0d0/sr3+acc*2.0d0)/auc) + 1

  do sc = 1, ns
     fsc = dble(sc) - 1.5d0
     tc = t*tab*fsc
     phic = phiab*fsc
     dc = d*(1.0d0+dab*fsc)
     xc = 0.5d0*dc*dcos(phic)
     yc = 0.5d0*dc*dsin(phic)
     zc = tc
     xyzc(1,sc) = xc
     xyzc(2,sc) = yc
     xyzc(3,sc) = zc

     rabc(1,sc) = phic
     rabc(2,sc) = tc/t

     do sn = 1, ns
        fsn = dble(sn) - 1.5d0
        un = 0
        do un1 = -unmax, unmax
           fun1 = dble(un1)
           do un2 = -unmax, unmax
              if(un1.ne.0.or.un2.ne.0.or.sn.ne.sc) then
                 fun2 = dble(un2)
                 phin = a1phi*fun1 + a2phi*fun2 + phiab*fsn
                 if(phin.gt.-pi+eps.and.phin.le.pi+eps) then
                    tn = t*(a1t*fun1+a2t*fun2+tab*fsn)
                    dn = d*(1.0d0+dab*fsn)
                    xn = 0.5d0*dn*dcos(phin)
                    yn = 0.5d0*dn*dsin(phin)
                    zn = tn
                    x = xn - xc
                    y = yn - yc
                    z = zn - zc
                    r = dsqrt(x**2+y**2+z**2)
                    if( r.le.cutlength ) then
                       un = un + 1
                       if(un.gt.nu) then
                          write(*,100) un
                          stop
                       endif

                       xyzn(1,un,sn,sc) = xn
                       xyzn(2,un,sn,sc) = yn
                       xyzn(3,un,sn,sc) = zn

                       rabn(1,un,sn,sc) = phin
                       rabn(2,un,sn,sc) = tn / t
                    endif
                 endif
              endif
           enddo
        enddo
        nabf(sn,sc) = un
     enddo
  enddo

100 format(" >> tbcoord >> out of range >> ",i8)
end subroutine tbcoord
