
!********************************************************************
!     calculate ec,ev,wx,wc:
!     ec    : renormilized conduction band
!     ev    : renormilized valence band
!     wc    : coulomb interaction
!     wx    : exchange interaction
!********************************************************************
subroutine iniga(muc,kc,nmkr,nmq,con,con11,con12,con21,&
     vq,vqc11,vqc12,vqc21,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,&
     ec,ev,ecs,evs,wc,wx,kcv,fcv)
  use common,only:nmrel,qmesh
  use common,only:nmukr,murmx,krmx
  use common,only:nc,nt,mumx,dt
  use common,only:bound
  use common,only:u0
  use common,only:wfqp,wfmesh
  use common,only:kcapex,kvapex
  use common,only:ns

  use exkataif,only : ftrans, qpencal, checkk, checkmu, checkrange

  implicit none
  !input
  integer muc,kc
  type(nmrel),INTENT(in)::nmkr
  type(qmesh)::nmq
  real*8,dimension(:,:),pointer::con,con11,con12,con21
  complex*16,dimension(:,:,:,:),pointer::vq,vqc11,vqc12,vqc21
  type(wfqp)::wfcoe1,wfcoe2,wfcoe3
  type(wfmesh)::nmwf2,nmwf3
  !output
  complex*16 ec(2,bound),ev(2,bound),&
       wc(2,2,bound,bound),wx(2,2,bound,bound)
  real*8 evs(2,bound),ecs(2,bound)
  integer kcv(bound,2,2,2)
  integer kcvpl(murmx,krmx,2,2,2)
  complex*16 fcv(bound,2,2,2)
  !local
  !number index
  integer nk1,nk2,kpoint
  !wave vector
  integer, dimension(2)::mur1,kr1,mucon1,kcon1,muval1,kval1,&
       mur2,kr2,mucon2,kcon2,muval2,kval2
  integer kcir,kt
  !ec and ev
  integer mucont,kcont,muvalt,kvalt
  complex*16 cct(2),cvt(2),qev(2,2),qec(2,2)
  real*8,dimension(2):: ec1,ev1
  complex*16,dimension(2)::ec2,ev2
  complex*16,dimension(2,ns,ns)::ec2arry,ev2arry
  integer muqin,qin
  real*8 muq,q
  complex*16 vqc(2,2)
  !wc and wx
  real*8 conc,cont
  complex*16 vqx(2,2)
  complex*16,dimension(2,2,ns,ns)::wcarry,wxarry
  !wave function coefficient (ec,ev,wc,wx)
  complex*16,dimension(2,2)::cc1,cv1,cc2,cv2
  !con and val wave vector check
  integer,dimension(2)::muconcom,muvalcom,kconcom,kvalcom
  logical cond
  !condition
  integer mubt2,muup2,kbt2,kup2,mubt3,muup3,kbt3,kup3
  !others
  integer i,j,i1,j1,i2,j2,m1,m2,sc,sn
  real*8 kvec(2)
  real*8 ui,vi
  character*2 type
  real*8,parameter::pi=3.1415926d0
  character*3 sub

  type='ga'
  sub='ini'

  if (nmukr.gt.bound) then
     write (*,100) bound
     stop
  endif

!!!!! FR for exchange CB !!!!!!!!!!!!!!
  muqin=muc
  qin=kc
!  muq=muqin*1.0d0
!  q=qin*1.0d0
  call ftrans(muqin,qin,2,vqx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!! QP energy !!!!!!!!!!!!!!!!!
  mubt2=nmwf2%mubt
  muup2=nmwf2%muup
  kbt2=nmwf2%kbt
  kup2=nmwf2%kup

  mubt3=nmwf3%mubt
  muup3=nmwf3%muup
  kbt3=nmwf3%kbt
  kup3=nmwf3%kup

  nk1=0

  nlineloop:do i=1,nmkr%nline

     mur1(:)=(/nmkr%mu(i),-nmkr%mu(i)/)
     mucon1=(mur1+muc)/2
     muval1=(mur1-muc)/2

     call checkmu(type,mucon1,muval1,i,nmkr)

     krloop:do j=1,nmkr%nk(i)

        nk1=nk1+1

        kr1(:)=(/nmkr%k(i,j),-nmkr%k(i,j)/)

        kcon1=(kr1+kc)/2
        kval1=(kr1-kc)/2

        kcv(nk1,1,1,:)=muval1(:)
        kcv(nk1,1,2,:)=mucon1(:)

        kcv(nk1,2,1,:)=kval1(:)
        kcv(nk1,2,2,:)=kcon1(:)

        kcvpl(i,j,:,:,:)=kcv(nk1,:,:,:)


        call checkk(type,kcon1,kval1,i,j,nmkr)

        do i1 = 1, 2
           select case(i1)
           case(1)
              kvapex = 'K1'
              kcapex = 'K1'
           case(2)
              kvapex = 'K2'
              kcapex = 'K2'
           end select

           kcir = muval1(i1)
           kt = kval1(i1)

           if (kvapex.eq.'K1') then

              call checkrange(kcir,mubt2,muup2,kt,kbt2,kup2,kvapex,sub)

              ev1(i1)=wfcoe2%ev(kcir,kt)
              cv1(i1,:)=wfcoe2%cv(kcir,kt,:)

              fcv(nk1,1,i1,:)=cv1(i1,:)

           else

              call checkrange(kcir,mubt3,muup3,kt,kbt3,kup3,kvapex,sub)

              ev1(i1)=wfcoe3%ev(kcir,kt)
              cv1(i1,:)=wfcoe3%cv(kcir,kt,:)

              fcv(nk1,1,i1,:)=cv1(i1,:)

           endif

           kcir=mucon1(i1)
           kt=kcon1(i1)

           if (kcapex.eq.'K1') then

              call checkrange(kcir,mubt2,muup2,kt,kbt2,kup2,kvapex,sub)

              ec1(i1)=wfcoe2%ec(kcir,kt)
              cc1(i1,:)=wfcoe2%cc(kcir,kt,:)

              fcv(nk1,2,i1,:)=cc1(i1,:)

           else

              call checkrange(kcir,mubt3,muup3,kt,kbt3,kup3,kvapex,sub)

              ec1(i1)=wfcoe3%ec(kcir,kt)
              cc1(i1,:)=wfcoe3%cc(kcir,kt,:)

              fcv(nk1,2,i1,:)=cc1(i1,:)

           endif

        enddo


        ec2arry(:,:,:)=dcmplx(0.0d0,0.0d0)
        ev2arry(:,:,:)=dcmplx(0.0d0,0.0d0)


        do i1=1,2
           select case(i1)
           case(1)
              kvapex='K1'
              kcapex='K1'
           case(2)
              kvapex='K2'
              kcapex='K2'
           end select
           mucont=mucon1(i1)
           kcont=kcon1(i1)
           muvalt=muval1(i1)
           kvalt=kval1(i1)
           cct(:)=cc1(i1,:)
           cvt(:)=cv1(i1,:)

           call  qpencal(mucont,kcont,muvalt,kvalt,cct,cvt,&
                vq,con,nmq,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,qev,qec)

           ev2arry(i1,:,:)=qev(:,:)
           ec2arry(i1,:,:)=qec(:,:)
        enddo

        ev2(:)=(/SUM(ev2arry(1,:,:)),SUM(ev2arry(2,:,:))/)
        ec2(:)=(/SUM(ec2arry(1,:,:)),SUM(ec2arry(2,:,:))/)

        ev2=ev2/(nc*nt)
        ec2=ec2/(nc*nt)

        ec(:,nk1)=ec1(:)+ec2(:)
        ev(:,nk1)=ev1(:)+ev2(:)

        !added
        ecs(:,nk1)=ec1(:)
        evs(:,nk1)=ev1(:)
        !

     enddo krloop
  enddo nlineloop

  nmukr=nk1

  !write (*,300)

  deallocate(vq)
  deallocate(con)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!Direct and Exchange CB!!!!!!!!!!!!

  nk1=0
  do i1=1,nmkr%nline

     mur1(:)=(/nmkr%mu(i1),-nmkr%mu(i1)/)

     mucon1=(mur1+muc)/2
     muval1=(mur1-muc)/2

     call checkmu(type,mucon1,muval1,i1,nmkr)

     do j1=1,nmkr%nk(i1)
        nk1=nk1+1

        kr1(:)=(/nmkr%k(i1,j1),-nmkr%k(i1,j1)/)

        kcon1=(kr1+kc)/2
        kval1=(kr1-kc)/2

        call checkk(type,kcon1,kval1,i1,j1,nmkr)

        do m1=1,2
           select case(m1)
           case(1)
              kvapex='K1'
              kcapex='K1'
           case(2)
              kvapex='K2'
              kcapex='K2'
           end select

           kcir=muval1(m1)
           kt=kval1(m1)

           if (kvapex.eq.'K1') then

              call checkrange(kcir,mubt2,muup2,kt,kbt2,kup2,kvapex,sub)

              cv1(m1,:)=wfcoe2%cv(kcir,kt,:)
           else

              call checkrange(kcir,mubt3,muup3,kt,kbt3,kup3,kvapex,sub)

              cv1(m1,:)=wfcoe3%cv(kcir,kt,:)
           endif

           kcir=mucon1(m1)
           kt=kcon1(m1)

           if (kcapex.eq.'K1') then

              call checkrange(kcir,mubt2,muup2,kt,kbt2,kup2,kvapex,sub)

              cc1(m1,:)=wfcoe2%cc(kcir,kt,:)
           else

              call checkrange(kcir,mubt3,muup3,kt,kbt3,kup3,kvapex,sub)

              cc1(m1,:)=wfcoe3%cc(kcir,kt,:)
           endif

        enddo


        nk2=0

        do i2=1,nmkr%nline

           mur2(:)=(/nmkr%mu(i2),-nmkr%mu(i2)/)

           mucon2=(mur2+muc)/2
           muval2=(mur2-muc)/2


           call checkmu(type,mucon2,muval2,i2,nmkr)

           do j2=1,nmkr%nk(i2)

              nk2=nk2+1

              kr2(:)=(/nmkr%k(i2,j2),-nmkr%k(i2,j2)/)

              kcon2=(kr2+kc)/2
              kval2=(kr2-kc)/2

              call checkk(type,kcon2,kval2,i2,j2,nmkr)

              do m1=1,2
                 select case(m1)
                 case(1)
                    kvapex='K1'
                    kcapex='K1'
                 case(2)
                    kvapex='K2'
                    kcapex='K2'
                 end select

                 kcir=muval2(m1)
                 kt=kval2(m1)
                 if (kvapex.eq.'K1') then

                    call checkrange(kcir,mubt2,muup2,kt,kbt2,kup2,kvapex,sub)

                    cv2(m1,:)=wfcoe2%cv(kcir,kt,:)
                 else

                    call checkrange(kcir,mubt3,muup3,kt,kbt3,kup3,kvapex,sub)

                    cv2(m1,:)=wfcoe3%cv(kcir,kt,:)
                 endif

                 kcir=mucon2(m1)
                 kt=kcon2(m1)
                 if (kcapex.eq.'K1') then

                    call checkrange(kcir,mubt2,muup2,kt,kbt2,kup2,kvapex,sub)

                    cc2(m1,:)=wfcoe2%cc(kcir,kt,:)
                 else

                    call checkrange(kcir,mubt3,muup3,kt,kbt3,kup3,kvapex,sub)

                    cc2(m1,:)=wfcoe3%cc(kcir,kt,:)
                 endif

              enddo

              wcarry=0.0d0
              wxarry=0.0d0

              do m1=1,2
                 do m2=1,2

                    muqin=mur2(m2)-mur1(m1)
                    qin=kr2(m2)-kr1(m1)

!!!!
                    if (mod(muqin,2).ne.0) then
                       write (*,400)
                       stop
                    endif

                    if (mod(qin,2).ne.0) then
                       write (*,400)
                       stop
                    endif
!!!


                    if (m1.eq.m2) then
                       vqc(:,:)=vqc11(:,:,muqin,qin)
                       conc=1.0d0+con11(muqin,qin)
                    endif

                    if (m1.eq.1.and.m2.eq.2) then
                       vqc(:,:)=vqc12(:,:,muqin,qin)
                       conc=1.0d0+con12(muqin,qin)
                    endif

                    if (m1.eq.2.and.m2.eq.1) then
                       vqc(:,:)=vqc21(:,:,muqin,qin)
                       conc=1.0d0+con21(muqin,qin)
                    endif

                    do sc=1,ns
                       do sn=1,ns

                          wcarry(m1,m2,sn,sc)=vqc(sn,sc)/conc*&
                               conjg(cc1(m1,sn))*cc2(m2,sn)*cv1(m1,sc)*conjg(cv2(m2,sc))

                          !if use 2K for FT
                          wxarry(m1,m2,sn,sc)=vqx(sn,sc)&
                          *conjg(cc1(m1,sc))*cv1(m1,sc)*cc2(m2,sn)*conjg(cv2(m2,sn))
                          !if use -2K for FT
                          !wxarry(m1,m2,sn,sc)=vqx(sn,sc)&
                          !     *conjg(cc1(m1,sn))*cv1(m1,sn)*cc2(m2,sc)*conjg(cv2(m2,sc))

                          ! NOTICE:
                          ! Exchange term of Eq.(15) in MGM a1176 is 
                          ! NOT correct. dielectric function is NOT needed.
                       enddo
                    enddo


                    wcarry=wcarry/(nc*nt)
                    wxarry=wxarry/(nc*nt)

                    wc(m1,m2,nk1,nk2)=SUM(wcarry(m1,m2,:,:))
                    wx(m1,m2,nk1,nk2)=SUM(wxarry(m1,m2,:,:))

                 enddo
              enddo

           enddo
        enddo
     enddo
  enddo


100 format (">>ini>>Hamiltonian matrix element exceeds",i5,&
       "l should be reduced" )
300 format("finish QP energy calculation")
400 format(">>iniga-dia>>wavevector for relative motion is wrong ")
500 format (">>iniga-dia>>wavevector range error")


  !deallocate(wfcoe1%ev,wfcoe2%ev,wfcoe3%ev)
  !deallocate(wfcoe1%ec,wfcoe2%ec,wfcoe3%ec)
  !deallocate(wfcoe1%cv,wfcoe2%cv,wfcoe3%cv)
  !deallocate(wfcoe1%cc,wfcoe2%cc,wfcoe3%cc)
  !deallocate(vqc11,vqc12,vqc21)
  !deallocate(con11,con12,con21)

  return
end subroutine iniga
