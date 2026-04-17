
subroutine ftdia(muc,kc,nmkr,nmk1,nmk2,nmq,&
     nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,con,con11,&
     con12,con21,vq,vqc11,vqc12,vqc21)

  use common,ONLY:nmrel,nmmesh,qmesh
  use common,only:murmx,krmx
  use common,only:ns,mu2d,mumx
  use common,only:u0,kapp
  use common,only:nc,nt
  use common,only:wfqp
  use common,only:wfmesh

  use exkataif,only : wflinecbi, wfcoecbi, mukrg, wflinecb, wfcoecb, ftrans, diacb, wfline, wfcoe, dia

  implicit none

  !input
  integer muc,kc
  type(nmrel),INTENT(in)::nmkr
  type(nmmesh)::nmk1,nmk2
  type(qmesh)::nmq
  type(wfmesh)::nmwf1,nmwf2,nmwf3
  type(wfqp)::wfcoe1,wfcoe2,wfcoe3,wfcoe4

  !output
  real*8,dimension(:,:),pointer::con,con11,con12,con21
  complex*16,dimension(:,:,:,:),pointer::vq,vqc11,vqc12,vqc21

  !local
  !ft and cons for direct CB
  integer nline1,nline2,nline,kpoint
  integer,dimension(mu2d)::mut,btt,upt,mu1,bt1,up1,mu2,bt2,up2
  integer,dimension(mumx)::mu,bt,up
  integer mul,mur,ktl,ktr
  integer,dimension(:),pointer::muarr,btarr,uparr
  !ft and cons for QB
  integer qpoint,mubt,muup,qbt,qup
  integer mucont,kcont,muvalt,kvalt
  real*8 con1
  !ft and cons for both CB and QP
  integer muqin,qin
  integer muq,q
  complex*16,dimension(2,2)::vq1,vq2
  !initial CB parameters
  !allocation
  integer svqc11,svqc12,svqc21,scon11,scon12,scon21
  logical cond
  !others
  integer i,j,sn,sc,ii

!!!!! FR and dielectric constants for direct CB !!!!!
  !write (*,100)
  nline1=nmkr%nline
  nline2=nmkr%nline
  do i=1,nline1
     kpoint=nmkr%nk(i)
     mut(i)=nmkr%mu(i)
     btt(i)=MINVAL(nmkr%k(i,1:kpoint))
     upt(i)=MAXVAL(nmkr%k(i,1:kpoint))
  enddo

  !define k range and then calculate energy and
  !coefficient in this range
  !wfcoe1: K->K, wfcoe2:K'->K'
  !intra CBI
  call wflinecbi(nmk1,nmk2,nmwf1,nmwf2)
  call wfcoecbi(nmwf1,nmwf2,wfcoe1,wfcoe2)

  do j=1,3
     select case (j)
     case(1)
        mu1=mut
        bt1=btt
        up1=upt
        mu2=mu1
        bt2=bt1
        up2=up1
     case(2)
        mu2=mut
        bt2=btt
        up2=upt
        mu1=-mu2
        bt1=-up2
        up1=-bt2
     case(3)
        mu1=mut
        bt1=btt
        up1=upt
        mu2=-mu1
        bt2=-up1
        up2=-bt1
     end select

     !estimate boundary for vqc and con
     call mukrg(nline1,mu1,bt1,up1,nline2,mu2,bt2,up2,nline,mu,bt,up)
     allocate (muarr(nline),btarr(nline),uparr(nline))

     do i=1,nline
        muarr(i)=mu(i)
        btarr(i)=bt(i)
        uparr(i)=up(i)
     enddo
     mul=MINVAL(muarr)
     mur=MAXVAL(muarr)
     ktl=MINVAL(btarr)
     ktr=MAXVAL(uparr)

     if (j.eq.1) then
        allocate(vqc11(2,2,mul:mur,ktl:ktr),stat=svqc11)
        allocate(con11(mul:mur,ktl:ktr),stat=scon11)
        cond=svqc11.eq.0.and.scon11.eq.0
     elseif (j.eq.2) then
        allocate(vqc12(2,2,mul:mur,ktl:ktr),stat=svqc12)
        allocate(con12(mul:mur,ktl:ktr),stat=scon12)
        cond=svqc12.eq.0.and.scon12.eq.0
     elseif (j.eq.3) then
        allocate(vqc21(2,2,mul:mur,ktl:ktr),stat=svqc21)
        allocate(con21(mul:mur,ktl:ktr),stat=scon21)
        cond=svqc21.eq.0.and.scon21.eq.0
     endif

     if (.not.cond) then
        write (*,500)
        stop
     endif

     nmwf1%mubt=idint(mul/2.0d0)-1
     nmwf1%muup=idint(mur/2.0d0)+1
     nmwf1%kbt=idint(ktl/2.0d0)-1
     nmwf1%kup=idint(ktr/2.0d0)+1
     nmwf1%nline=nmwf1%muup-nmwf1%mubt+1
     nmwf1%point=nmwf1%kup-nmwf1%kbt+1

     !define k range and then calculate energy and
     !coefficient in this range
     !wfcoe3: K->K', wfcoe4:K'->K
     !For interCBI
     call wflinecb(nmk1,nmk2,nmwf1,nmwf2,nmwf3)
     call wfcoecb(nmwf1,nmwf2,nmwf3,wfcoe3,wfcoe4)

     do i=1,nline
        muqin=mu(i)
        do qin=bt(i),up(i),1
           if (mod(muqin,2).ne.0) cycle
           if (mod(qin,2).ne.0) cycle

           muq=idint(muqin/2.0d0)
           q=idint(qin/2.0d0)

           !1:not include on-site
           !2:include on-site
           call ftrans(muq,q,2,vq2)

           vq2=vq2/kapp

           call diacb(muq,q,nmk1,nmk2,nmwf2,nmwf3,&
                wfcoe1,wfcoe2,wfcoe3,wfcoe4,vq2,con1)

           if (j.eq.1) then
              con11(muqin,qin)=con1
           elseif (j.eq.2) then
              con12(muqin,qin)=con1
           elseif (j.eq.3) then
              con21(muqin,qin)=con1
           endif

           do sc=1,ns
              do sn=1,ns
                 if (j.eq.1) then
                    vqc11(sn,sc,muqin,qin)=vq2(sn,sc)
                 elseif(j.eq.2) then
                    vqc12(sn,sc,muqin,qin)=vq2(sn,sc)
                 elseif(j.eq.3) then
                    vqc21(sn,sc,muqin,qin)=vq2(sn,sc)
                 endif
              enddo
           enddo
        enddo
     enddo

     deallocate(muarr,btarr,uparr)


     deallocate(wfcoe3%ev,wfcoe4%ev)
     deallocate(wfcoe3%ec,wfcoe4%ec)
     deallocate(wfcoe3%cv,wfcoe4%cv)
     deallocate(wfcoe3%cc,wfcoe4%cc)


  enddo

  deallocate(wfcoe1%ev,wfcoe2%ev)
  deallocate(wfcoe1%ec,wfcoe2%ec)
  deallocate(wfcoe1%cv,wfcoe2%cv)
  deallocate(wfcoe1%cc,wfcoe2%cc)

  !write (*,200)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!FT and dielectric constants for QP!!!
  !write (*,300)


  nline1=nmq%nline

  allocate (vq(1:2,1:2,nmq%mu(1):nmq%mu(nline1),nmq%bt:nmq%up))
  allocate (con(nmq%mu(1):nmq%mu(nline1),nmq%bt:nmq%up))
  nmwf1%nline=nmq%nline
  nmwf1%mubt=nmq%mu(1)
  nmwf1%muup=nmq%mu(nline1)
  nmwf1%kbt=nmq%bt
  nmwf1%kup=nmq%up
  nmwf1%point=nmwf1%kup-nmwf1%kbt+1


  !define k range and then calculate energy and
  !coefficient in this range
  !wfcoe1:Gamma (qline), wfcoe2:Gamma+K, wfcoe3:Gamma+K'
  !For QP
  call wfline(nmk1,nmk2,nmwf1,nmwf2,nmwf3)
  call wfcoe(nmwf1,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3)

  lineloop:do i=1,nmq%nline
     muqin=nmq%mu(i)
     qloop:do qin=nmq%bt,nmq%up,1

        call ftrans(muqin,qin,1,vq1)

        do sc=1,ns
           do sn=1,ns
              vq(sn,sc,muqin,qin)=vq1(sn,sc)/kapp
           enddo
        enddo

        vq2=vq1
        vq2(1,1)=vq2(1,1)+u0
        vq2(2,2)=vq2(2,2)+u0

        vq2=vq2/kapp

        call dia(muqin,qin,nmk1,nmk2,nmwf2,nmwf3,&
             wfcoe1,wfcoe2,wfcoe3,vq2,con1)

        con(muqin,qin)=con1

     enddo qloop
  enddo lineloop
  !write (*,400)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

100 format("begin FT and dielectric constants for direct CB")
200 format("end FT and dielectric constants for direct CB")
300 format("begin FT and dielectric constants for QP")
400 format("end FT and dielectric constants for QP")
500 format (">>ftdia.f90>>allocation error ")

  return
end subroutine ftdia
