
!= exkatampi.f90
! * Excitation energy kataura plot based on ETB
!
!=== Notice
! * module common must be complied first
!
!=== Input
! * strrel.dat : parameter file for ETB based on tbdftsnt.
! * exkataparameter.f90 : parameter for exkatmpi.f90
! * see INPUT PARAMETER in this source code
!
!=== Output
! * the output file name is defined by variable fn3.
! 1. n
! 2. m
! 3. diameter [nm]
! 4. chiral angle [rad]
! 5. excitation energy [eV]
! 6. binding energy [eV]
! 7. quasi particle energy [eV]
! 8. self energy [eV]
! 9. self energy - binding energy [eV]
!
!=== References
! 1. J. Jiang <em>et al</em>., Phys. Rev. B 75, 035407 (2007).
! Chirality dependence of the exciton effects in Single-Wall Carbon Nanotubes: Tight-binding model
!
!=== Contact
! author:: Jie Jiang
! date::  Apr 26, 2006
! e-mail:: mailto:jiang@chips.ncsu.edu
!
! modify and comment:: Kentaro Sato
! e-mail:: mailto:kentaro@flex.phys.tohoku.ac.jp

module common
  integer,parameter :: mu2d=50, mumx=102, murmx=50, krmx=1000, bound=2000
  integer nc,nt,nmukr
  real(8) kbase1(2),kbase2(2)
  real(8) v0

  real(8) :: u0, kapp
  ! u0 and kapp: see INPUT PARAMETER in this source code

  type::nmmesh
     integer nline
     integer,dimension(mu2d)::mu,bt,up,vhsk
     real(8),dimension(mu2d)::vhse,vhsec,vhsev
  end type nmmesh

  type::qmesh
     integer nline
     integer,dimension(mumx)::mu
     integer bt,up,point
  end type qmesh

  type::wfmesh
     integer nline,mubt,muup
     integer kbt,kup,point
  end type wfmesh

  type::wfqp
     real(8),pointer::ev(:,:),ec(:,:)
     complex(8),pointer:: cv(:,:,:),cc(:,:,:)
  end type wfqp

  type::nmcen
     integer nline
     integer,dimension(mumx)::mu,bt,up
  end type nmcen

  type::nmrel
     integer nline
     integer,dimension(murmx)::mu,muc,muv,nk
     integer,dimension(murmx,krmx)::k,kc,kv
  end type nmrel

  integer,parameter::nu=64000,ns=2,nu1=1000
  integer nabf(ns,ns)
  real(8) xyzc(3,ns)
  real(8) xyzn(3,nu,ns,ns)
  real(8) rabc(2,ns)
  real(8) rabn(2,nu,ns,ns)
  !parameter for tube
  real(8) tt,dt,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi
  integer nab(ns,ns)
  real(8) rab(3,64,ns,ns)
  real(8) hab(4,4,64,ns,ns)
  real(8) oab(4,4,64,ns,ns)

  integer unmax
  real(8) rabc1(2,ns)
  real(8) rabn1(2,-nu1:nu1,-nu1:nu1,ns,ns)
  !apex for 2D BZ
  character*2 kcapex,kvapex

  integer nn,mm

  complex(8) :: kexchange(4), kdirect(4)

end module common



!------------------------------------
! include subroutines and functions
!------------------------------------

include 'exkataparameter.f90'
include 'exkataif.f90'
include 'eg3.f90'
include 'dsort.f90'
include 'num2str.f90'
include 'arguments.f90'
include 'numvhs.f90'

include 'tubpar.f'
include 'cutline.f90'
include 'tbtube.f'
include 'tbtubecb.f'
include 'tbpot.f'
include 'wfline.f90'
include 'wfcoe.f90'
include 'findkc.f90'
include 'findkr.f90'
include 'iniga-dia.f90'
include 'mukrg.f90'
include 'ftrans.f90'
include 'ftranscb.f90'
include 'tbcoord.f90'
include 'exhs.f90'
include 'dia.f90'
include 'diacb.f90'
include 'qpencal.f90'
include 'cutoff.f90'
include 'ftdia.f90'
include 'wflinecbi.f90'
include 'wfcoecbi.f90'
include 'wflinecb.f90'
include 'wfcoecb.f90'
include 'check.f90'
include 'qline.f90'
include 'inik.f90'


program exkatampi
  != Excitation energy kataura plot based on ETB

  use common

  use exkataif,only : cutline, qline, ftdia, iniga, exh, findkc, findkr, num2str, numvhs, inik
  use exkataparameter,only : PI, AUC

  implicit none

  include 'mpif.h'

  ! input. see INPUT PARAMETER in this source code.
  character(len=1) :: stype
  character(len=2) :: otype, tubetype
  character(len=3) :: ao
  integer :: cline, vline, kc, state, subband
  real(8) :: length

  ! output
  integer nd
  integer kcv(bound,2,2,2)
  real(8)  ecvc,ecvsc
  complex(8) fcv(bound,2,2,2)
  type(nmcen)::nmkc

  ! local

  !tube parameters
  real(8) kcx,kcy,ktx,kty
  real(8) dtinnm,dtina
  real(8) vab(64,ns,ns)

  !type
  type(nmmesh)::nmk1,nmk2
  type(nmrel)::nmkr
  type(wfqp)::wfcoe1,wfcoe2,wfcoe3
  type(wfmesh)::nmwf2,nmwf3
  type(qmesh)::nmq
  !for ftdia and iniga
  real(8),dimension(:,:),pointer::con,con11,con12,con21
  complex(8),dimension(:,:,:,:),pointer::vq,vqc11,vqc12,vqc21
  !h elements
  complex(8) ec(2,bound),ev(2,bound),&
       wc(2,2,bound,bound),wx(2,2,bound,bound)
  real(8) ecs(2,bound),evs(2,bound)
  !findkc and findkr
  character*2 type
  !h type and results
  real(8) ex(bound)
  complex(8) exvec(bound,bound)
  !other
  integer i,j
  real(8),allocatable:: ect(:),evt(:),ecst(:),evst(:)
  real(8) ecc,evc,ecsc,evsc
  complex(8) vec(bound)

  integer muc

  integer :: tubetypeindex
  real(8) :: dtmin, dtmax
  ! for tube selection. parameter(dtmin=0.46,dtmax=1.6)

  integer :: iterat
  real(8) :: lini, en, toler
  ! for reading ETB parameter file named ``fn10''
  
  real(8) :: ttini, dtini, tabini, phiabini, dabini
  ! for subroutine tubpar

  character(len=100) :: fn3, fn10, fn20
  ! fn3  : output file name
  ! fn10 : input file name
  ! see INPUT PARAMETER in this source code

  integer(4) :: ierr, p, my_rank, rank, dest, source, count, countq
  integer,parameter :: tag = 1, tagq = 2, tagke = 3, tagkd = 4, tagw1 = 5, tagw2 = 6, tagw22 = 7, tagc1 = 8, tagc2 = 9, tagc22 = 10, tagpos = 11, tagwcn2g = 12
  integer(4) :: status(MPI_STATUS_SIZE)
  integer(4),parameter :: size = 7, sizeq = 2, sizeke = 4, sizekd = 4, sizew1 = 1, sizew2 = 1, sizec1 = 1, sizec2 = 1, sizepos = 4, sizewcn2g = 1
  real(8) :: mxs(size), mxr(size), qvecs(2), qvecr(2)
  ! for MPI

  integer :: io
  ! file status

  integer,parameter :: AsciiZero = 48
  ! Ascii number 48 correspond to 0. Reference. man ascii.

  integer :: b

  integer :: muq, kq, lq
  character(4) :: senbu

  integer :: MaxNumVHS

  ! for short-range Coulomb parameter w1 and w2
  complex(8) :: wc2(2,2,bound,bound), wx2(2,2,bound,bound), wc3(2,2,bound,bound), wx3(2,2,bound,bound), srcp(4)
  integer :: pos1, pos2

  complex(8) :: wcn2g


  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  call SRArguments(my_rank,p,tubetype,subband)

  !--------------------------
  ! START : INPUT PARAMETER
  !--------------------------

  u0 = 11.3d0
  ! the energy cost to place two electrons on a single site
  ! see eq.(14), J. Jiang et al, Phys. Rev. B 75, 035407 (2007).

  kapp = 2.22d0
  ! a static dielectric constant describing the effects of electrons in core states
  ! see eq.(6) and the following paragraph, J. Jiang et al, Phys. Rev. B 75, 035407 (2007).

  !tubetype = 's0'
  if(tubetype.eq.'s1') tubetypeindex = 1
  if(tubetype.eq.'s2') tubetypeindex = 2
  if(tubetype.eq.'s0') tubetypeindex = 0
  ! nanotube type
  ! s1 : Semiconductor Type I
  ! s2 : Semiconductor Type II
  ! s0 : Metal

  dtmin = 0.7D0
  ! dtmin : minimum diameter [nm]
  dtmax = 2.0D0
  ! dtmax : maximum diameter [nm]
  ! diameter range depends on information of nanotube structure optimazation by ETB.
  ! this information is stored in "strrel.dat".
  ! if you want to calculate large diameter nanotube, you run relaxte.f90 and get information.

  otype = 'go'
  ! orbital type
  ! go : A2
  ! ge : A1
  ! k1 : E
  ! k2 : E*

  stype = 's'
  ! sping type
  ! s : singlet
  ! t : triplet

  state = 1
  ! exciton level
  ! n = 1, 2,...
  ! WARNING: n=1 correspond to exciton level 0.

  length = 200.0d0
  ! nanotube length [nm]

  ao = 'one'
  ! all : all cutting lines
  ! one : one cutting lines

  cline = subband
  ! electron cutting line 

  vline = subband
  ! hole cutting line

  kc = 0
  ! 1D center of mass momentum

!  fn3 = './def/qex-kat-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
  fn3 = './results/qex-kat-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
       &//char(state-1+AsciiZero)//'-'//char(idint(kapp*1.0d0)+48)//'.dat'
  ! output filename

!  fn10 = '/home1/students/kentaro/for/datafiles/etb/tbdftsnt/strrel.dat'
  fn10 = './strrel.dat'
  ! input data from ETB

  !------------------------
  ! END : INPUT PARAMETER
  !------------------------


  rank = 0
  if( my_rank.ne.0 ) then
     ! calculate exciton energy kataura plot for one nanotube and
     ! send results to rank 0 machine

     open(21,file=fn10,status='old',iostat=io)
     if (io.ne.0) then
        write(*,401)fn10
        stop
     endif

     do i=1,5
        read(21,100,iostat=io)
        if (io.ne.0) then
           write(*,402)fn10
           stop
        endif
     enddo

     do
        read(21,115,iostat=io) nn,mm,lini,tt,dt,theta,tab,phiab,dab,en,toler,iterat
        if (io.ne.0) exit

        if (mod(2*nn+mm,3).ne.tubetypeindex) cycle

        if( (otype=='go').or.(otype=='ge') ) then
           MaxNumVHS = numvhs(nn,mm)
           if( MaxNumVHS < subband ) cycle
        end if
        
        dtinnm=dt
        dtina=dt/auc

        if (dtinnm.lt.dtmin.or.dtinnm.gt.dtmax) cycle

        !initial data for ETB wavefunction
        call tubpar (nn,mm,length,ttini,dtini,theta,tabini,phiabini,&
             dabini,a1t,a1phi,a2t,a2phi,nc,nt,kcx,kcy,ktx,kty)
        call tbtini (tt,dt,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nab,&
             rab,hab,oab,vab)

        rank=rank+1
        if(rank.ge.p) rank=1
        if (my_rank.eq.rank) then

           kbase1=(/kcx,kcy/)
           kbase2=(/ktx,kty/)

           !tube coordinate
           call tbcoord
           !cutting lines around K
           call cutline(1,nmk1)
           !cutting lines around K'
           call cutline(2,nmk2)
           !cutting lines around G
           call qline(nmk1,nmq)

           muc = nmk1%mu(cline) - nmk1%mu(vline)

           select case(otype)
           case('ge')
              type = 'ga'
           case('go')
              type = 'ga'
           case('k1')
              type = 'k1'
           case('k2')
              type = 'k2'
           end select

           !find center-of-mass momentum kc-kv
           call findkc(type,nmk1,nmk2,ao,cline,vline,nmkc)
           !find relative motion momentum kc+kv
           call findkr(type,nmk1,nmk2,muc,kc,ao,cline,vline,nmkr)

           !dielectric function (con, con**)
           !con:QP con11,con12,con21:CB
           !11:K-->K or K'-->K', 12:K-->K', 21:K'-->K
           !CB potential FT (vq, vq**)
           !vq:QP vq11,vq12,vq21:CB
           !11:K-->K or K'-->K', 12:K-->K', 21:K'-->K
           call ftdia(muc,kc,nmkr,nmk1,nmk2,nmq,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,con,con11,con12,con21,vq,vqc11,vqc12,vqc21)

           ! write dielectric function
           senbu =  num2str(nn,mm)
           !  fn20 = './results/def-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
           fn20 = './def/def-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
                &//char(state-1+AsciiZero)//'-'//char(idint(kapp*1.0d0)+48)//'-swnt'//senbu//'.dat'
           open(20,FILE=fn20)
           do lq = 1, nmq%nline
              muq = nmq%mu(lq)
              if( muq /= 0 ) cycle
              do kq = nmq%bt, nmq%up
                 write(20,*) kq*((2.0d0*pi/tt)/nt), 1.0D0+con(muq,kq)
              end do
              write(20,*)
           end do
           close(20)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(type.eq.'ga') then
              !A excitons
              !ec:conduction QP, ev:valence QP, ecs:conduction SP
              !ecv:valence SP
              !wc:Kd, wx:Kx
              call iniga(muc,kc,nmkr,nmq,con,con11,con12,con21,vq,vqc11,vqc12,vqc21,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,ec,ev,ecs,evs,wc,wx,kcv,fcv,wc2,wx2,wc3,wx3)
           else if( (type=='k1').or.(type=='k2') ) then
              call inik(muc,kc,nmkr,ec,ev,wc,wx)
           end if

           pos1 = (nmk1%vhsk(cline)-1)*nmk2%nline + 1
           pos2 = pos1 + nmk2%vhsk(vline) - 1
           srcp(1) = wx(1,1,pos1,pos2)
           srcp(2) = wx(2,1,pos1,pos2)
           srcp(3) = wx(1,2,pos1,pos2)
           srcp(4) = wx(2,2,pos1,pos2)

           nd=nmukr

           !exciton eigenvalues and eigenvectors
           !ex(state):velues exvec(:,state):vectors
           call exh(otype,stype,nd,nmkr,ec,ev,wc,wx,wc2,wx2,wc3,wx3,wcn2g,ex,exvec)

           vec(:)=exvec(:,state)

           allocate (ect(nd),evt(nd))
           allocate (ecst(nd),evst(nd))
           do j=1,nd
              ect(j)=real(ec(1,j))
              evt(j)=real(ev(1,j))
              ecst(j)=ecs(1,j)
              evst(j)=evs(1,j)
           enddo
           ecc=MINVAL(ect)
           evc=MAXVAL(evt)
           ecvc=ecc-evc
           deallocate (ect,evt)

           ecsc=MINVAL(ecst)
           evsc=MAXVAL(evst)
           ecvsc=ecsc-evsc
           deallocate (ecst,evst)

           !dt
           mxs(1)=dtinnm
           !theta
           mxs(2)=theta
           !excitation enegy
           mxs(3)=ex(1)
           !binding energy
           mxs(4)=ecvc-ex(1)
           !QP energy
           mxs(5)=ecvc
           !selfenergy
           mxs(6)=ecvc-ecvsc
           !selfenergy-bd
           mxs(7)=mxs(6)-mxs(4)
           ! q vector
           qvecs(:) = nmk1%mu(subband)*kbase1(:) + nmk1%vhsk(subband)*kbase2(:)

           dest=0
           call MPI_SEND (mxs,size,MPI_COMPLEX16,dest,tag,MPI_COMM_WORLD,ierr)
           call MPI_SEND (qvecs,sizeq,MPI_REAL8,dest,tagq,MPI_COMM_WORLD,ierr)
           call MPI_SEND (kexchange,sizeke,MPI_COMPLEX16,dest,tagke,MPI_COMM_WORLD,ierr)
           call MPI_SEND (kdirect,sizekd,MPI_COMPLEX16,dest,tagkd,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wx3(1,1,1,1),sizew1,MPI_COMPLEX16,dest,tagw1,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wx3(1,2,1,1),sizew2,MPI_COMPLEX16,dest,tagw2,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wx3(2,1,1,1),sizew2,MPI_COMPLEX16,dest,tagw22,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wc3(1,1,1,1),sizec1,MPI_COMPLEX16,dest,tagc1,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wc3(1,2,1,1),sizec2,MPI_COMPLEX16,dest,tagc2,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wc3(2,1,1,1),sizec2,MPI_COMPLEX16,dest,tagc22,MPI_COMM_WORLD,ierr)
           call MPI_SEND (srcp,sizepos,MPI_COMPLEX16,dest,tagpos,MPI_COMM_WORLD,ierr)
           call MPI_SEND (wcn2g,sizewcn2g,MPI_COMPLEX16,dest,tagwcn2g,MPI_COMM_WORLD,ierr)
        end if

     end do

     close (21,status='keep')

  end if


  if( my_rank.eq.0 ) then
     ! gather results from other machines and write results to output file

     open(53,file=fn3,status='unknown')
     close(53,status='delete')
     open (53,file=fn3,status='new')

     rank=0

     open(21,file=fn10,status='old',iostat=io)
     if (io.ne.0) then
        write(*,401)fn10
        stop
     endif

     do i=1,5
        read(21,100,iostat=io)
        if (io.ne.0) then
           write(*,402)fn10
           stop
        endif
     enddo

     do
        read(21,115,iostat=io) nn,mm,lini,tt,dt,theta,tab,phiab,dab,en,toler,iterat
        if (io.ne.0) exit

        if (mod(2*nn+mm,3).ne.tubetypeindex) cycle

        if( (otype=='go').or.(otype=='ge') ) then
           MaxNumVHS = numvhs(nn,mm)
           if( MaxNumVHS < subband ) cycle
        end if
        
        dtinnm = dt
        dtina = dt / auc
        if(dtinnm.lt.dtmin.or.dtinnm.gt.dtmax) cycle

        rank = rank + 1
        if(rank.ge.p) rank = 1
        source = rank
        call MPI_RECV (mxr,size,MPI_COMPLEX16,source,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (qvecr,sizeq,MPI_REAL8,source,tagq,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_REAL8, countq, ierr )
        call MPI_RECV (kexchange,sizeke,MPI_COMPLEX16,source,tagke,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (kdirect,sizekd,MPI_COMPLEX16,source,tagkd,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wx3(1,1,1,1),sizew1,MPI_COMPLEX16,source,tagw1,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wx3(1,2,1,1),sizew2,MPI_COMPLEX16,source,tagw2,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wx3(2,1,1,1),sizew2,MPI_COMPLEX16,source,tagw22,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wc3(1,1,1,1),sizec1,MPI_COMPLEX16,source,tagc1,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wc3(1,2,1,1),sizec2,MPI_COMPLEX16,source,tagc2,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wc3(2,1,1,1),sizec2,MPI_COMPLEX16,source,tagc22,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (srcp,sizepos,MPI_COMPLEX16,source,tagpos,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (wcn2g,sizewcn2g,MPI_COMPLEX16,source,tagwcn2g,MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )


        write(53,110) nn, mm, mxr(:), qvecr(:)&
             &, wcn2g, abs(wcn2g)&
             &, srcp(1:4), abs(srcp(1:4))&
             &, wx3(1,1,1,1), abs(wx3(1,1,1,1)), wx3(2,1,1,1), abs(wx3(2,1,1,1)), wx3(1,2,1,1), abs(wx3(1,2,1,1)), wx3(2,1,1,1)+wx3(1,2,1,1), abs(wx3(1,1,1,1)+wx3(1,2,1,1))&
             &, wc3(1,1,1,1), abs(wc3(1,1,1,1)), wc3(2,1,1,1), abs(wc3(2,1,1,1)), wc3(1,2,1,1), abs(wc3(1,2,1,1)), wc3(2,1,1,1)+wc3(1,2,1,1), abs(wc3(1,1,1,1)+wc3(1,2,1,1))&
             &, kexchange, abs(kexchange), 0.5D0*(kexchange(1)+kexchange(4)), abs(0.5D0*(kexchange(1)+kexchange(4)))&
             &, 0.5D0*(kexchange(2)+kexchange(3)), abs(0.5D0*(kexchange(2)+kexchange(3)))&
             &, kdirect, abs(kdirect), 0.5D0*(kdirect(1)+kdirect(4)), abs(0.5D0*(kdirect(1)+kdirect(4)))&
             &, 0.5D0*(kdirect(2)+kdirect(3)), abs(0.5D0*(kdirect(2)+kdirect(3)))
        ! write calculating results to the file named ``fn3''
        ! see INPUT in the head of this source code.

        write(*,FMT='(A10,2I4)') "finish: ", nn, mm
     end do

     close (21,status='keep')

  end if


110 format(2i5,9f12.5,75F15.10)
120 format(6f12.5)
115 format(" ",ss,i3.3," ",ss,i3.3," ",sp,d22.15," ",sp,d22.15," ",sp,d2&
       &2.15," ",sp,d22.15," ",sp,d22.15," ",sp,d22.15," ",sp,d22.15," ",s&
       &p,d22.15," ",sp,d22.15," ",ss,i6.6," ")
100 format(a223)
401 format(" file ",a10," is not found >> ")
402 format(" file ",a10," has wrong format >> &
       &")
502 format(" (",ss,i3.3,",",ss,i3.3,") ")

  call MPI_FINALIZE (ierr)

end program exkatampi
