
!*************************************************************************
!     calculate ec,ev,wx,wc:
!     ec    : renormilized conduction band
!     ev    : renormilized valence band
!     wc    : coulomb interaction
!     wx    : exchange interaction
!*************************************************************************

subroutine inik(muc,kc,nmkr,ec,ev,wc,wx)
!subroutine inik(muc,kc,nmkr,ec,ev,wc,wx,kcv,fcv)

  use exkataparameter,only : Egresult, PI
  use common,only:kbase1,kbase2
  use common,ONLY:nmrel,nmukr
  use common,only:nc,nt,ns
  use common,only:bound
  use exkataif,only : eg3
  use common,only : tt, nt

  implicit none

  !     input
  integer muc,kc
  type(nmrel),INTENT(in)::nmkr
  character*2 type
  !     output
  complex*16 ec(2,bound),ev(2,bound),&
       wc(2,2,bound,bound),wx(2,2,bound,bound)
  real*8 kcv(bound,2,2,2)
  complex*16 fcv(bound,2,2,2)
  !     local
  type (EGresult) ::e_3
  !     number index
  integer nk1,nk2
  !     wave vector
  integer, dimension(2)::mur1,kr1,mucon1,kcon1,muval1,kval1,&
       mur2,kr2,mucon2,kcon2,muval2,kval2
  !     ec and ev
  real*8,dimension(2):: ec1,ev1
  complex*16,dimension(2)::ec2,ev2
  complex*16,dimension(2,ns,ns)::ec2arry,ev2arry
  !    wc and wx
  complex*16,dimension(2,2,ns,ns)::wcarry,wxarry
  !    wave function coefficient (ec,ev,wc,wx)
  complex*16,dimension(2,2)::cc1,cv1,cc2,cv2
  !    fourier transformation
  integer muqin,qin
  real*8 muq,q
  complex*16,dimension(:,:,:,:),allocatable::vq
  complex*16,dimension(2,2)::vq1,vqe,vqc,vqx
  !   ft for ex
  integer nline
  integer,dimension(:),allocatable::mu,bt,up
  complex*16,dimension(:,:,:,:),allocatable::vqc1
  integer mul,mur,ktl,ktr
  !    con and val wave vector check
  integer,dimension(2)::muconcom,muvalcom,kconcom,kvalcom
  logical cond
  !    others
  integer i,j,i1,j1,i2,j2,m1,m2,sc,sn
  real*8 kvec(2)
  real(8) :: stepsize

  if (nmukr.gt.bound) then
     write (*,100) bound
     stop
  endif

  stepsize = (2.0d0*PI/tt)/nt

  allocate (vq(1:2,1:2,-nc/2:nc/2-1,-nt/2:nt/2-1))

  do muqin=-nc/2,nc/2-1
     do qin=-nt/2,nt/2-1

        muq=muqin*1.0d0
        q=qin*1.0d0

        call ftrans(muq,q,1,vq1)

        do sc=1,ns
           do sn=1,ns
              vq(sn,sc,muqin,qin)=vq1(sn,sc)
           enddo
        enddo

     enddo
  enddo


  muqin=2*muc
  qin=2*kc

  muq=muqin*1.0d0
  q=qin*1.0d0

  call ftrans(muq,q,2,vqx)

  nline=nmkr%nline
  allocate (mu(nline),bt(nline),up(nline))

  do i=1,nline
     mu(i)=nmkr%mu(i)
     bt(i)=MINVAL(nmkr%k(i,:))
     up(i)=MAXVAL(nmkr%k(i,:))
  enddo

  mul=MINVAL(mu)
  mur=MAXVAL(mu)
  ktl=MINVAL(bt)
  ktr=MAXVAL(up)

  allocate(vqc1(2,2,mul:mur,ktl:ktr))


  do i=1,nline
     do qin=bt(i),up(i),1

        muqin=mu(i)

        muq=muqin*0.5d0
        q=qin*0.5d0

        call ftrans(muq,q,2,vq1)

        do sc=1,ns
           do sn=1,ns
              vqc1(sn,sc,muqin,qin)=vq1(sn,sc)
           enddo
        enddo

     enddo
  enddo


  deallocate(mu,bt,up)


  nk1=0
  do i=1,nmkr%nline

     mur1(1)=nmkr%mu(i)

     mucon1(1)=(mur1(1)+muc)/2
     muval1(1)=(mur1(1)-muc)/2
     !
     muconcom(1)=nmkr%muc(i)
     muvalcom(1)=nmkr%muv(i)

     cond=mucon1(1).ne.muconcom(1).or.muval1(1).ne.muvalcom(1)
     if (cond) then
        write (*,50)
        stop
     endif
     !

     do j=1,nmkr%nk(i)

        nk1=nk1+1

        ec(:,nk1)=0.0d0
        ev(:,nk1)=0.0d0

        kr1(1)=nmkr%k(i,j)

        kcon1(1)=(kr1(1)+kc)/2
        kval1(1)=(kr1(1)-kc)/2
        !
        kconcom(1)=nmkr%kc(i,j)
        kvalcom(1)=nmkr%kv(i,j)
        cond=kcon1(1).ne.kconcom(1).or.kval1(1).ne.kvalcom(1)
        if (cond) then
           write (*,50)
           stop
        endif

        kcv(nk1,1,1,1:2)=(/mucon1(1)*1.0d0,kcon1(1)*stepsize/)
        kcv(nk1,2,1,1:2)=(/muval1(1)*1.0d0,kval1(1)*stepsize/)

        !
        do i1=1,1
           kvec=mucon1(i1)*kbase1+kcon1(i1)*stepsize*&
                kbase2/sqrt(dot_product(kbase2,kbase2))
           e_3=eg3(kvec(1),kvec(2))
           ec1(i1)=e_3%w2
           cc1(i1,1:2)=(/e_3%ca(2),e_3%cb(2)/)

           fcv(nk1,1,i1,1:2)=cc1(i1,1:2)

           kvec=muval1(i1)*kbase1+kval1(i1)*stepsize*&
                kbase2/sqrt(dot_product(kbase2,kbase2))
           e_3=eg3(kvec(1),kvec(2))
           ev1(i1)=e_3%w1
           cv1(i1,1:2)=(/e_3%ca(1),e_3%cb(1)/)

           fcv(nk1,2,i1,1:2)=conjg(cv1(i1,1:2))

        enddo

        ec2arry(:,:,:)=dcmplx(0.0d0,0.0d0)
        ev2arry(:,:,:)=dcmplx(0.0d0,0.0d0)

        do muqin=-nc/2,nc/2-1
           do qin=-nt/2,nt/2-1

              mucon2(1)=mucon1(1)+muqin
              kcon2(1)=kcon1(1)+qin
              muval2(1)=muval1(1)-muqin
              kval2(1)=kval1(1)-qin

              do i1=1,1
                 kvec=mucon2(i1)*kbase1+kcon2(i1)*stepsize*&
                      kbase2/sqrt(dot_product(kbase2,kbase2))
                 e_3=eg3(kvec(1),kvec(2))
                 cc2(i1,1:2)=(/e_3%ca(2),e_3%cb(2)/)
                 kvec=muval2(i1)*kbase1+kval2(i1)*stepsize*&
                      kbase2/sqrt(dot_product(kbase2,kbase2))
                 e_3=eg3(kvec(1),kvec(2))
                 cv2(i1,1:2)=(/e_3%ca(1),e_3%cb(1)/)
              enddo


              do sc=1,ns
                 do sn=1,ns
                    vqe(sn,sc)=vq(sn,sc,muqin,qin)
                 enddo
              enddo


              do i1=1,1
                 do sc=1,ns
                    do sn=1,ns

                       ec2arry(i1,sn,sc)=ec2arry(i1,sn,sc)+vqe(sn,sc)*conjg(cc1(i1,sn))*cc1(i1,sc)*&
                            cc2(i1,sn)*conjg(cc2(i1,sc))

                       ev2arry(i1,sn,sc)=ev2arry(i1,sn,sc)-vqe(sn,sc)*cv1(i1,sn)*conjg(cv1(i1,sc))*&
                            conjg(cv2(i1,sn))*cv2(i1,sc)

                    enddo
                 enddo
              enddo

           enddo
        enddo


        ec2arry=ec2arry/(nc*nt)
        ev2arry=ev2arry/(nc*nt)

        ec2(1)=SUM(ec2arry(1,:,:))
        ev2(1)=SUM(ev2arry(1,:,:))

        ec(1,nk1)=ec1(1)+ec2(1)
        ev(1,nk1)=ev1(1)+ev2(1)


     enddo
  enddo

  nmukr=nk1

  write (*,*) 'ok2'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  nk1=0
  do i1=1,nmkr%nline


     mur1(1)=nmkr%mu(i1)

     mucon1(1)=(mur1(1)+muc)/2
     muval1(1)=(mur1(1)-muc)/2
     !
     muconcom(1)=nmkr%muc(i1)
     muvalcom(1)=nmkr%muv(i1)
     cond=mucon1(1).ne.muconcom(1).or.muval1(1).ne.muvalcom(1)
     if (cond) then
        write (*,50)
        stop
     endif
     !

     do j1=1,nmkr%nk(i1)

        nk1=nk1+1

        kr1(1)=nmkr%k(i1,j1)

        kcon1(1)=(kr1(1)+kc)/2
        kval1(1)=(kr1(1)-kc)/2
        !
        kconcom(1)=nmkr%kc(i1,j1)
        kvalcom(1)=nmkr%kv(i1,j1)
        cond=kcon1(1).ne.kconcom(1).or.kval1(1).ne.kvalcom(1)
        if (cond) then
           write (*,50)
           stop
        endif
        !

        do m1=1,1

           kvec=mucon1(m1)*kbase1+kcon1(m1)*stepsize*&
                kbase2/sqrt(dot_product(kbase2,kbase2))

           e_3=eg3(kvec(1),kvec(2))
           cc1(m1,1:2)=(/e_3%ca(2),e_3%cb(2)/)

           kvec=muval1(m1)*kbase1+kval1(m1)*stepsize*&
                kbase2/sqrt(dot_product(kbase2,kbase2))

           e_3=eg3(kvec(1),kvec(2))
           cv1(m1,1:2)=(/e_3%ca(1),e_3%cb(1)/)

        enddo


        nk2=0
        do i2=1,nmkr%nline

           mur2(1)=nmkr%mu(i2)

           mucon2(1)=(mur2(1)+muc)/2
           muval2(1)=(mur2(1)-muc)/2
           !
           muconcom(1)=nmkr%muc(i2)
           muvalcom(1)=nmkr%muv(i2)
           cond=mucon2(1).ne.muconcom(1).or.muval2(1).ne.muvalcom(1)
           if (cond) then
              write (*,50)
              stop
           endif
           !

           do j2=1,nmkr%nk(i2)
              nk2=nk2+1

              wc(:,:,nk1,nk2)=0.0d0
              wx(:,:,nk1,nk2)=0.0d0

              kr2(1)=nmkr%k(i2,j2)

              kcon2(1)=(kr2(1)+kc)/2
              kval2(1)=(kr2(1)-kc)/2
              !
              kconcom(1)=nmkr%kc(i2,j2)
              kvalcom(1)=nmkr%kv(i2,j2)
              cond=kcon2(1).ne.kconcom(1).or.kval2(1).ne.kvalcom(1)
              if (cond) then
                 write (*,50)
                 stop
              endif
              !

              do m1=1,1

                 kvec=mucon2(m1)*kbase1+kcon2(m1)*stepsize*&
                      kbase2/sqrt(dot_product(kbase2,kbase2))

                 e_3=eg3(kvec(1),kvec(2))
                 cc2(m1,1:2)=(/e_3%ca(2),e_3%cb(2)/)

                 kvec=muval2(m1)*kbase1+kval2(m1)*stepsize*&
                      kbase2/sqrt(dot_product(kbase2,kbase2))

                 e_3=eg3(kvec(1),kvec(2))
                 cv2(m1,1:2)=(/e_3%ca(1),e_3%cb(1)/)

              enddo


              wcarry=0.0d0
              wxarry=0.0d0


              do m1=1,1
                 do m2=1,1


                    muqin=mur2(m2)-mur1(m1)
                    qin=kr2(m2)-kr1(m1)

                    vqc(:,:)=vqc1(:,:,muqin,qin)

                    do sc=1,ns
                       do sn=1,ns

                          wcarry(m1,m2,sn,sc)=&
                               vqc(sn,sc)*conjg(cc1(m1,sn))*cc2(m2,sn)*cv1(m1,sc)*conjg(cv2(m2,sc))

                          wxarry(m1,m2,sn,sc)=&
                               vqx(sn,sc)*conjg(cc1(m1,sc))*cv1(m1,sc)*cc2(m2,sn)*conjg(cv2(m2,sn))

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


50 format (">>ini>>con or val wave vector is not correct")
100 format (">>ini>>Hamiltonian matrix element exceeds",i5,"l should be reduced")

  deallocate(vq,vqc1)

end subroutine inik
