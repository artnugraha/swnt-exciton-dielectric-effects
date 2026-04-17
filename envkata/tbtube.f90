!>
!! \file tbtube.f90
!! \brief tight-binding dispersion for the rolled up carbon nanotube
!<

!>
!! \defgroup tbdftsnt Extended Tight Binding module
!<

!>
!! \ingroup tbdftsnt
!! \brief tight-binding dispersion for the rolled up carbon nanotube
!!
!! \note
!! - matrix blocks
!!  - bn=1:nb -- (2s,2px,2py,2pz) atomic orbital at neighbor site
!!  - bc=1:nb -- (2s,2px,2py,2pz) atomic orbital at central site
!!  - dt=1:nd -- (t,phi,d) cylindrical coordinates
!!  - un=1:nab(sn,sc) -- neighbor unit cell index
!!  - sn=1:ns -- (A,B) type of neighbor site
!!  - sc=1:ns -- (A,B) type of central site
!!
!! \param[in] t  nanotube unit cell length [nm]
!! \param[in] d  nanotube diameter [nm]
!! \param[in] tab  translational distance between A & B atoms [t]
!! \param[in] phiab  circumferential distance between A & B atoms [radian]
!! \param[in] dab  radial distance between A & B atoms [d]
!! \param[in] a1t  translational component of graphene unit vector a1 [t]
!! \param[in] a1phi  circumferential component of a1 vector [radian]
!! \param[in] a2t  translational component of graphene unit vector a2 [t]
!! \param[in] a2phi  circumferential component of a2 vector [radian]
!!
!! \param[out] nab(ns,ns) number of atom pairs of a given type
!! \param[out] rab(nd,nu,ns,ns) interatomic distances for each atom pair
!! \param[out] hab(nb,nb,nu,ns,ns) hamiltonian matrix blocks
!! \param[out] oab(nb,nb,nu,ns,ns) overlap matrix blocks
!! \param[out] vab(nu,ns,ns) repulsive potentials between atom sites
!!
!! \section References
!! -# V. N. Popov, New J. Phys. 6, 17 (2004).
!!  - http://dx.doi.org/10.1088/1367-2630/6/1/017
!! -# Ge. G. Samsonidze et al., Appl. Phys. Lett. 85, 5703 (2004).
!!  - http://dx.doi.org/10.1063/1.1829160
!!
!! \author Created by Georgii Samsonidze <gsm@mgm.mit.edu>
!! \author Commented by Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!!
!! \date last modified : 2004-07-02
!! \date date created : 2004-05-11
!<
subroutine tbtini(t,d,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nab,rab,hab,oab,vab)

  implicit none

  integer,parameter :: no = 3, nb = 4, nd = 3, nu = 64, ns = 2
  real(8),parameter :: pi = +0.3141592653589793d+01
  real(8),parameter :: sr3 = +0.1732050807568877d+01
  real(8),parameter :: acc = +0.142d+00
  real(8),parameter :: auc = acc*sr3
  real(8),parameter :: bohr = +0.5291772108d-01
  real(8),parameter :: cutoff = +0.1d+02*bohr
  real(8),parameter :: eps = +0.1d-08

  ! Input
  real(8),intent(in) :: t, d, tab, phiab, dab, a1t, a1phi, a2t, a2phi

  ! Output
  integer,intent(out) :: nab(ns,ns)
  real(8),intent(out) :: rab(nd,nu,ns,ns)
  real(8),intent(out) :: hab(nb,nb,nu,ns,ns)
  real(8),intent(out) :: oab(nb,nb,nu,ns,ns)
  real(8),intent(out) :: vab(nu,ns,ns)

  ! Local
  integer :: sc,sn,un,unmax,un1,un2,i,j,k,l
  real(8) :: cylcut,fsc,fsn,fun1,fun2
  real(8) :: tc,phic,dc,tn,phin,dn,xc,yc,zc,xn,yn,zn,x,y,z,r
  real(8) :: xr,yr,zr,xti,yti,zti,xto,yto,zto,lr,lti,lto
  real(8) :: pcn(nd,nd),pco(nd,nd),pno(nd,nd)
  real(8) :: pci(nd),pnf(nd)
  real(8) :: tbt(no,no),tbs(no,no)
  real(8) :: tbv

  if (d.gt.cutoff) then
     cylcut=d*dasin(cutoff/d)
  elseif (d.lt.cutoff) then
     cylcut=dsqrt((d*pi/2.0d0)**2+(cutoff**2-d**2))
  else
     cylcut=d*pi/2.0d0
  endif

  unmax=idint((cylcut*2.0d0/sr3+acc*2.0d0)/auc)+1
  do sc=1,ns
     fsc=dfloat(sc)-1.5d0
     tc=t*tab*fsc
     phic=phiab*fsc
     dc=d*(1.0d0+dab*fsc)
     xc=0.5d0*dc*dcos(phic)
     yc=0.5d0*dc*dsin(phic)
     zc=tc
     do sn=1,ns
        fsn=dfloat(sn)-1.5d0
        un=0
        do un1=-unmax,unmax
           fun1=dfloat(un1)
           do un2=-unmax,unmax
              if (un1.ne.0.or.un2.ne.0.or.sn.ne.sc) then
                 fun2=dfloat(un2)
                 phin=a1phi*fun1+a2phi*fun2+phiab*fsn
                 if (phin.gt.-pi+eps.and.phin.le.pi+eps) then
                    tn=t*(a1t*fun1+a2t*fun2+tab*fsn)
                    dn=d*(1.0d0+dab*fsn)
                    xn=0.5d0*dn*dcos(phin)
                    yn=0.5d0*dn*dsin(phin)
                    zn=tn
                    x=xn-xc
                    y=yn-yc
                    z=zn-zc
                    r=dsqrt(x**2+y**2+z**2)
                    if (r.le.cutoff) then
                       un=un+1
                       if (un.gt.nu) then
                          write(*,100)
                          stop
                       endif
                       rab(1,un,sn,sc)=(tn-tc)/t
                       rab(2,un,sn,sc)=(phin-phic)
                       rab(3,un,sn,sc)=(dn-dc)/d
                       call tbpot (r,tbt,tbs,tbv)
                       vab(un,sn,sc)=tbv

                       xr=xn-xc
                       yr=yn-yc
                       zr=zn-zc
                       lr=dsqrt(xr**2+yr**2+zr**2)
                       if (lr.lt.eps) then
                          write(*,111)
                          stop
                       endif
                       xr=xr/lr
                       yr=yr/lr
                       zr=zr/lr
                       if (dabs(xr).gt.eps.or.dabs(yr).gt.eps) then
                          if (yr.gt.+0.0d+00) then
                             xto=+yr
                             yto=-xr
                             zto=+0.0d+00
                          else
                             xto=-yr
                             yto=+xr
                             zto=+0.0d+00
                          endif
                       else
                          xto=xn+xc
                          yto=yn+yc
                          zto=+0.0d+00
                          lto=dsqrt(xto**2+yto**2+zto**2)
                          if (lto.lt.eps) then
                             xto=+0.1d+01
                             yto=+0.0d+00
                             zto=+0.0d+00
                          endif
                       endif
                       lto=dsqrt(xto**2+yto**2+zto**2)
                       if (lto.lt.eps) then
                          write(*,113)
                          stop
                       endif
                       xto=xto/lto
                       yto=yto/lto
                       zto=zto/lto
                       xti=yto*zr-zto*yr
                       yti=zto*xr-xto*zr
                       zti=xto*yr-yto*xr
                       lti=dsqrt(xti**2+yti**2+zti**2)
                       if (lti.lt.1.0d0-eps.or.lti.gt.1.0d0+eps) then
                          write(*,112)
                          stop
                       endif
                       xti=xti/lti
                       yti=yti/lti
                       zti=zti/lti

                       pcn(1,1)=xr
                       pcn(2,1)=yr
                       pcn(3,1)=zr
                       pcn(1,2)=xti
                       pcn(2,2)=yti
                       pcn(3,2)=zti
                       pcn(1,3)=xto
                       pcn(2,3)=yto
                       pcn(3,3)=zto

                       pco(1,1)=dcos(phic)
                       pco(2,1)=dsin(phic)
                       pco(3,1)=0.0d0
                       pco(1,2)=-dsin(phic)
                       pco(2,2)=dcos(phic)
                       pco(3,2)=0.0d0
                       pco(1,3)=0.0d0
                       pco(2,3)=0.0d0
                       pco(3,3)=1.0d0

                       pno(1,1)=dcos(phin)
                       pno(2,1)=dsin(phin)
                       pno(3,1)=0.0d0
                       pno(1,2)=-dsin(phin)
                       pno(2,2)=dcos(phin)
                       pno(3,2)=0.0d0
                       pno(1,3)=0.0d0
                       pno(2,3)=0.0d0
                       pno(3,3)=1.0d0

                       hab(1,1,un,sn,sc)=+tbt(1,1)
                       oab(1,1,un,sn,sc)=+tbs(1,1)

                       do j=2,nb
                          pnf(1)=0.0d0
                          do l=1,nd
                             pnf(1)=pnf(1)+pcn(l,1)*pno(l,j-1)
                          enddo
                          hab(j,1,un,sn,sc)=+tbt(2,1)*pnf(1)
                          oab(j,1,un,sn,sc)=+tbs(2,1)*pnf(1)
                       enddo

                       do i=2,nb
                          pci(1)=0.0d0
                          do l=1,nd
                             pci(1)=pci(1)+pcn(l,1)*pco(l,i-1)
                          enddo
                          hab(1,i,un,sn,sc)=+tbt(1,2)*pci(1)
                          oab(1,i,un,sn,sc)=+tbs(1,2)*pci(1)
                       enddo

                       do i=2,nb
                          do j=2,nb
                             do k=1,nd
                                pci(k)=0.0d0
                                pnf(k)=0.0d0
                             enddo
                             do k=1,nd
                                do l=1,nd
                                   pci(k)=pci(k)+pcn(l,k)*pco(l,i-1)
                                   pnf(k)=pnf(k)+pcn(l,k)*pno(l,j-1)
                                enddo
                             enddo
                             hab(j,i,un,sn,sc)=+tbt(2,2)*pci(1)*pnf(1)+&
                                  &tbt(3,3)*pci(2)*pnf(2)+&
                                  &tbt(3,3)*pci(3)*pnf(3)
                             oab(j,i,un,sn,sc)=+tbs(2,2)*pci(1)*pnf(1)+&
                                  &tbs(3,3)*pci(2)*pnf(2)+&
                                  &tbs(3,3)*pci(3)*pnf(3)
                          enddo
                       enddo

                    endif
                 endif
              endif
           enddo
        enddo
        nab(sn,sc)=un
     enddo
  enddo

100 format(" >> tbdftsnt >> tbtube >> tbtini >> nab exceeds nu >> ")
111 format(" >> tbdftsnt >> tbtube >> tbtini >> lr >> out of range >> ")
112 format(" >> tbdftsnt >> tbtube >> tbtini >> lti >> out of range >> ")
113 format(" >> tbdftsnt >> tbtube >> tbtini >> lto >> out of range >> ")

end subroutine tbtini



!>
!! \ingroup tbdftsnt
!! \brief tight-binding dispersion for the rolled up carbon nanotube
!!
!! \note
!! - Eigen values is ascending order for the specification of LAPACK. 
!! - tight binding parameters have been changed by Georgii-san
!!  - OLD: e2s=-13.573[eV] & e2p=-4.882[eV] & ef=-4.24778676[eV]
!!  - NEW: e2s=-13.573[eV] & e2p=-5.372[eV] & ef=-4.7504448646575[eV]
!! - matrix blocks
!!  - bn=1:nb -- (2s,2px,2py,2pz) atomic orbital at neighbor site
!!  - bc=1:nb -- (2s,2px,2py,2pz) atomic orbital at central site
!!  - dt=1:nd -- (t,phi,d) cylindrical coordinates
!!  - un=1:nab(sn,sc) -- neighbor unit cell index
!!  - sn=1:ns -- (A,B) type of neighbor site
!!  - sc=1:ns -- (A,B) type of central site
!!
!! \param[in] evflag option for LAPACK zhegv
!! - evflag is 0 : calculate eigen values
!! - evflag is not 0 : calculate eigen values and eigen vectors
!! \param[in] kc wavevector in circumferential direction (1-nc/2...nc/2)
!! \param[in] kt wavevector in translational direction (1-nt/2...nt/2)
!! \param[in] nc number of cutting lines in circumferential direction
!! \param[in] nt number of mesh points in translational direction
!! \param[in] nab number of atom pairs [dimensionless]
!! \param[in] rab interatomic distances [t,radian,d]
!! \param[in] hab hamiltonian matrix blocks [eV]
!! \param[in] oab overlap matrix blocks [eV]
!!
!! \param[out] ener eigenvalue (or band energy) [eV]
!! - ener[4] is pi band around K or K'point
!! - ener[5] is anti-pi band around K or K'point
!! \param[out] ampl wave function coefficient
!! - ampl(i,j) is i-th component of normalized j-th eigenvector
!!
!! \section References
!! -# V. N. Popov, New J. Phys. 6, 17 (2004).
!!  - http://dx.doi.org/10.1088/1367-2630/6/1/017
!! -# Ge. G. Samsonidze et al., Appl. Phys. Lett. 85, 5703 (2004).
!!  - http://dx.doi.org/10.1063/1.1829160
!!
!! \author Created by Georgii Samsonidze <gsm@mgm.mit.edu>
!! \author Commented by Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!!
!! \date last modified : 2004-07-02
!! \date date created : 2004-05-11
!<
subroutine tbtube(evflag,kc,kt,nc,nt,nab,rab,hab,oab,ener,ampl)

  implicit none

  integer,parameter :: ne = 2, nb = 4, nd = 3, nu = 64, ns = 2
  integer,parameter :: nm = nb*ns

!  real(8),parameter :: onsite(ne) = (/ -0.13573d+02, -0.4882d+01 /)
!  real(8),parameter :: ef = -0.424778676d+01
  real(8),parameter :: onsite(ne) = (/ -0.13573d+02, -0.5372d+01 /)
  real(8),parameter :: ef = -0.475044486465750d+01

  real(8),parameter :: pi = +0.3141592653589793d+01
  complex(8),parameter :: c0 = ( +0.0d+00, +0.0d+00 )
  complex(8),parameter :: c1 = ( +0.1d+01, +0.0d+00 )

  ! Input
  integer,intent(in) :: evflag, kc, kt, nc, nt
  integer,intent(in) :: nab(ns,ns)
  real(8),intent(in) :: rab(nd,nu,ns,ns)
  real(8),intent(in) :: hab(nb,nb,nu,ns,ns)
  real(8),intent(in) :: oab(nb,nb,nu,ns,ns)

  ! Output
  real(8),intent(out) :: ener(nm)
  complex(8),intent(out) :: ampl(nm,nm)

  ! Local
  integer :: i, j, bc, bn, sc, sn, un
  real(8) :: fkt, fkphi
  complex(8) :: pf
  complex(8) :: h(nm,nm), o(nm,nm)

  !=begin=lapack=data=!
  integer,parameter :: n = nm
  character(len=1),parameter :: ceval = 'N', cevec = 'V', uplo = 'U'
  integer,parameter :: itype = 1, lda = n, ldb = n, lwork = (32+1)*n, lrwork = 3*n-2
  character(len=1) :: jobz
  integer :: info
  complex(8) :: a(lda,n)
  complex(8) :: b(ldb,n)
  real(8) :: w(n)
  complex(8) :: work(lwork)
  real(8) :: rwork(lrwork)
  !=end=lapack=data=!

  fkt = dfloat(kt)/dfloat(nt)*2.0d0*pi
  fkphi = dfloat(kc)

  do i=1,ns
     h(1+(i-1)*nb,1+(i-1)*nb)=dcmplx(onsite(1),0.0d0)
     do j=2,nb
        h(j+(i-1)*nb,j+(i-1)*nb)=dcmplx(onsite(2),0.0d0)
     enddo
  enddo

  do i=1,nm
     do j=1,i-1
        h(j,i)=c0
     enddo
     do j=i+1,nm
        h(j,i)=c0
     enddo
  enddo

  do i=1,nm
     o(i,i)=c1
  enddo

  do i=1,nm
     do j=1,i-1
        o(j,i)=c0
     enddo
     do j=i+1,nm
        o(j,i)=c0
     enddo
  enddo

  do sc=1,ns
     bc=nb*(sc-1)
     do sn=1,ns
        bn=nb*(sn-1)
        do un=1,nab(sn,sc)
           !pf=cdexp(dcmplx(0.0d0,-fkt*rab(1,un,sn,sc)-fkphi*rab(2,un,sn,sc)))
           pf=cdexp(dcmplx(0.0d0,fkt*rab(1,un,sn,sc)+fkphi*rab(2,un,sn,sc)))

           do i=1,nb
              do j=1,nb
                 h(j+bn,i+bc)=h(j+bn,i+bc)+pf*dcmplx(hab(j,i,un,sn,sc),0.0d0)
              enddo
           enddo
           do i=1,nb
              do j=1,nb
                 o(j+bn,i+bc)=o(j+bn,i+bc)+pf*dcmplx(oab(j,i,un,sn,sc),0.0d0)
              enddo
           enddo
        enddo
     enddo
  enddo


  !=begin=lapack=code=!
  if (evflag.eq.0) then
     jobz = ceval
  else
     jobz = cevec
  end if

  a(1:n,1:n) = h(1:n,1:n)
  b(1:n,1:n) = o(1:n,1:n)

  call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
  if( info.ne.0 ) then
     ! diagonalization is faild
     write(*,100)
     stop
  end if

  ener(1:n) = w(1:n) - ef

  if( evflag /= 0 ) then
     ampl(1:n,1:n) = a(1:n,1:n)
  end if
  !=end=lapack=code=!


100 format(" >> tbdftsnt >> tbtube >> lapack >> zhegv >> info >> overflow >> ")

end subroutine tbtube
