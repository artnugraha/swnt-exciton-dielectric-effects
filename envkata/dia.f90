!>
!! \file dia.f90
!! \brief calculate Coulomb potential
!<

!>
!! \brief calculate multiplying Coulomb potentail by polarization
!!
!! \param[in] x is energy (conduction or valence band)
!!
!! \param[out] con multiplying Coulomb potentail by polarization
!!
!! \section References
!! -# Eq.(7), J. Jiang et al., Phys.Rev.B 75, 035407 (2007).
!!  - http://link.aps.org/abstract/PRB/v75/e035407
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
subroutine dia(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,vq2,con)

  use common,only : nc, nt, nmmesh, wfqp, wfmesh
  use exkataif,only : self

  implicit none

  ! input
  integer :: muq, q
  complex(8) :: vq2(2,2)
  type(nmmesh) :: nmk1, nmk2
  type(wfmesh) :: nmwf2, nmwf3
  type(wfqp) :: wfcoe1, wfcoe2, wfcoe3

  ! output
  real(8),intent(out) :: con

  ! local
  complex(8) :: vq(2,2), vqt
  real(8) :: selfenergy
  integer :: sn, sc

  vqt = 0.0d0
  do sn = 1, 2
     do sc = 1, 2
        vqt = vqt + vq2(sn,sc)
     end do
  end do
  vqt = vqt / 4.0d0
  ! Eq.(18), J. Jiang et al., Phys. Rev. B 75, 035407 (2007).

  call self(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,selfenergy)
  con = real(vqt)*selfenergy/(nc*nt)
  ! multiplying Coulomb potentail by polarization in Eq.(7)

end subroutine dia



subroutine self(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,selfenergy)
  != calculate polarization in Eq.(17), J. Jiang et al., Phys. Rev. B 75, 035407 (2007).
  !
  !=== Input
  !=== Output

  use common,only:nc,nt
  use common,only:nmmesh
  use common,only:wfmesh
  use common,only:wfqp

  use exkataif,only : cutoff !, checkrange

  !use common,only:nab,rab,hab,oab
  implicit none

  !input
  integer muq,q
  type(nmmesh)::nmk1,nmk2
  type(wfmesh)::nmwf2,nmwf3
  type(wfqp)::wfcoe1,wfcoe2,wfcoe3

  !output
  real(8) selfenergy

  !local
  !wavevector
  integer muk,k
  integer kc1,kc2,kt1,kt2
  !ec and ev
  real(8) ev1,ec1,ev2,ec2
  !wave function coefficient
  complex(8),dimension(2)::cv1,cc1,cv2,cc2
  !tbtube
  integer*4 nb,nd,nu,ns
  parameter (nb=4,nd=3,nu=64,ns=2)
  integer*4 nm
  parameter (nm=nb*ns)
  integer*4 evflag
  parameter (evflag=1)
  real(8) ener(nm)
  complex(8) ampl(nm,nm)
  !pol
  !integer s,o1,o2
  real(8) f1d,f2d
  complex(8) f1u,f2u
  real(8) gv1,gv2,gc1,gc2,g1,g2,pol
  !cutoff
!  real(8),external::cutoff
  !for check
  integer mubt2,muup2,mubt3,muup3,kbt2,kup2,kbt3,kup3
  logical cond1,cond2
  !
  integer i

  mubt2=nmwf2%mubt
  muup2=nmwf2%muup
  kbt2=nmwf2%kbt
  kup2=nmwf2%kup

  mubt3=nmwf3%mubt
  muup3=nmwf3%muup
  kbt3=nmwf3%kbt
  kup3=nmwf3%kup

  pol=0.0d0

  do i = 1, nmk1%nline
     muk = nmk1%mu(i)
     do k = nmk1%bt(i), nmk1%up(i), 1
        kc1 = muk
        kt1 = k

        kc2 = muk + muq
        kt2 = k + q

        call checkrange(kc1,mubt2,muup2,kt1,kbt2,kup2)
        call checkrange(kc2,mubt2,muup2,kt2,kbt2,kup2)

        ev1=wfcoe2%ev(kc1,kt1)
        ec1=wfcoe2%ec(kc1,kt1)
        cv1(:)=wfcoe2%cv(kc1,kt1,:)
        cc1(:)=wfcoe2%cc(kc1,kt1,:)

        ev2=wfcoe2%ev(kc2,kt2)
        ec2=wfcoe2%ec(kc2,kt2)
        cv2(:)=wfcoe2%cv(kc2,kt2,:)
        cc2(:)=wfcoe2%cc(kc2,kt2,:)

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



  do i = 1, nmk2%nline
     muk = nmk2%mu(i)
     do k = nmk2%bt(i), nmk2%up(i), 1
        kc1 = muk
        kt1 = k

        kc2 = muk + muq
        kt2 = k + q


        call checkrange(kc1,mubt3,muup3,kt1,kbt3,kup3)
        call checkrange(kc2,mubt3,muup3,kt2,kbt3,kup3)

        ev1=wfcoe3%ev(kc1,kt1)
        ec1=wfcoe3%ec(kc1,kt1)
        cv1(:)=wfcoe3%cv(kc1,kt1,:)
        cc1(:)=wfcoe3%cc(kc1,kt1,:)

        ev2=wfcoe3%ev(kc2,kt2)
        ec2=wfcoe3%ec(kc2,kt2)
        cv2(:)=wfcoe3%cv(kc2,kt2,:)
        cc2(:)=wfcoe3%cc(kc2,kt2,:)


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

100 format (">>dia>>wavevector range is not correct for K1")
200 format (">>dia>>wavevector range is not correct for K2")

end subroutine self
