
!>
!! \file envkata.f90
!>

!>
!! \mainpage
!!
!! \section Input
!! -# Input Parameters
!!  - variables for exkatampi system.

!! -# strrel.dat
!!  - parameter file for ETB. 
!!    you must use this, otherwise you cannot get the family behaviour. 
!!    Please see Georgii-san's PhD thesis.
!! -# exkataparameter.f90
!!  - parameters for exkatampi system.
!!    Probably you need not change parameters in exkataparameter.f90.
!!
!! \section Output
!! Output file is kapp-sX-eXX-XXX.dat and qexkat-sX-eXX-XXX.dat
!! - Filename
!!  - sX is 
!!   - Metal : s0
!!   - Type I : s1
!!   - Type II : s2
!!  - eXX is Eii
!!  - XXX is kappa value: 010 = 1.0, 020 = 2.0, etc
!!   - A1 : ge
!!   - A2 : go
!! - Meaning of each column in qexkat-sX-eXX-XXX.dat
!!  -# n
!!  -# m
!!  -# diameter [nm]
!!  -# chiral angle [rad]
!!  -# exciton energy [eV] \f$ \Omega _n \f$
!!  -# binding energy [eV] \f$ E_{\rm bd} = E_{\rm QP} - \Omega _0 \f$
!!  -# quasi particle energy [eV] \f$ E({\bf k}_{\rm c}) =
!!  \epsilon({\bf k}_{\rm c}) + \Sigma({\bf k}_{\rm c}) \f$
!!  -# self energy [eV] \f$ \Sigma = \sum W \f$
!!  -# self energy - binding energy [eV] \f$ E_{\rm mb} =
!!  \Sigma - E_{\rm bd} \f$
!!  -# lk [invers nm]
!!  -# polarization function at lk
!!  -# kappa
!!  -# qvecx
!!  -# qvecy
!!
!! \attention
!! -# module common must be compiled first
!! -# This program is based on MGM a1176. However some equations are 
!!   incorrect. Please be careful. KS wrote comment in this source code or 
!!   in PhD thesis.
!! -# Original programs have been developped using Intel Fortran 7.1.
!!    Probably you need IFC 7.1.
!!
!! \section subroutine and function
!! - \ref tbdftsnt
!!
!! \section References
!! -# J. Jiang et al., Phys. Rev. B 75, 035407 (2007)
!!  - http://link.aps.org/abstract/PRB/v75/e035407
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!! \author A.R.T. Nugraha (nugraha@flex.phys.tohoku.ac.jp)
!<

!>
!! \brief main program for Exciton Energy Kataura Plot
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!! \author A.R.T. Nugraha (nugraha@flex.phys.tohoku.ac.jp)
!<
program main

  use common
  use exkataif,only: cutline, qline, ftdia, iniga, exh, &
       & findkc, findkr, num2str, numvhs, inik
  use exkataparameter,only: PI, AUC

  implicit none

  include 'mpif.h'

  !
  ! Input Parameters. Check Input Parameters in exkatampi.f90.
  !

  !> spin type
  !! - s : singlet
  !! - t : triplet
  character(len=1) :: stype

  !> orbital type. select exciton type.
  !! - go : A2
  !! - ge : A1
  !! - k1 : E
  !! - k2 : E*
  character(len=2) :: otype

  !> nanotube type
  !! - s1 : Semiconductor Type I
  !! - s2 : Semiconductor Type II
  !! - s0 : Metal
  character(2) :: tubetype

  !> How many cutting lines do you take into calculation?
  !! Generally one cutting line is enough.
  !! - one : one cutting lines
  !! - all : all cutting lines
  character(3) :: ao

  !> electron cutting line
  !! - E11 : 1
  !! - E22 : 2
  !! - Eii : i
  integer :: cline

  !> hole cutting line
  !! - E11 : 1
  !! - E22 : 2
  !! - Eii : i
  integer :: vline

  !> nanotube length. default is 200nm. 
  real(8) :: length

  !> 1D center of mass momentum
  integer :: kc

  !> exciton level
  !! WARNING: state = 1 corresponds to excitation level 0.
  !! - state = 1: excitation level is 0
  !! - state = 2: excitation level is 1
  !! - state = i: excitation level is i-1
  integer :: state

  !> maximum diameter [nm]
  real(8) :: dtmin
  !> minimum diameter [nm]
  real(8) :: dtmax

  !> fn3  : output file name (qex-...)
  character(100) :: fn3
  !> fn4
  character(100) :: fn4
  !! fn10 : input parameter file name (strrel.dat)
  character(100) :: fn10
  character(100) :: fn5
  ! fn3  : output file name
  ! fn4  : output file name
  ! fn5  : output file name
  ! fn10 : input file name
  ! see INPUT PARAMETER in this source code

  integer :: subband

  !
  ! Output
  !
  real(8) :: ecvc, ecvsc
  real(8) :: lk, polareps

  !
  ! Local
  !
  integer :: tubetypeindex
  type(nmmesh) :: nmk1, nmk2
  type(nmrel) :: nmkr
  type(wfqp) :: wfcoe1, wfcoe2, wfcoe3
  type(wfmesh) :: nmwf2, nmwf3
  type(nmcen) :: nmkc
  integer :: kcv(bound,2,2,2)
  integer :: nd
  complex(8) :: fcv(bound,2,2,2)


  !tube parameters
  real(8) :: kcx, kcy, ktx, kty
  real(8) :: dtinnm, dtina
  real(8) :: vab(64,ns,ns)

  !type
  type(qmesh) :: nmq
  !for ftdia and iniga
  real(8),dimension(:,:),pointer :: con, con11, con12, con21
  complex(8),dimension(:,:,:,:),pointer :: vq, vqc11, vqc12, vqc21
  !h elements
  complex(8) :: ec(2,bound),ev(2,bound), wc(2,2,bound,bound), &
       & wx(2,2,bound,bound)
  real(8) :: ecs(2,bound), evs(2,bound)
  !findkc and findkr
  character*2 :: type
  !h type and results
  real(8) :: ex(bound)
  complex(8) :: exvec(bound,bound)
  !other
  integer :: i, j
  real(8),allocatable :: ect(:), evt(:), ecst(:), evst(:)
  real(8) :: ecc, evc, ecsc, evsc
  complex(8) vec(bound)

  integer :: muc

  integer :: iterat
  real(8) :: lini, en, toler
  ! for reading ETB parameter file named ``fn10''
  
  real(8) :: ttini, dtini, tabini, phiabini, dabini
  ! for subroutine tubpar

  integer(4) :: ierr, p, my_rank, rank, dest, source, count, countq, countk
  integer,parameter :: tag = 1, tagq = 2, tagk = 3
  integer(4) :: status(MPI_STATUS_SIZE)
  integer(4),parameter :: size = 9, sizeq = 2, sizek = 1
  real(8) :: mxs(size), mxr(size), qvecs(2), qvecr(2)
  ! for MPI

  integer :: io
  ! file status

  ! Ascii number 48 correspond to character 0.
  ! Try to "man ascii" on your console.
  integer,parameter :: AsciiZero = 48

  integer :: b

  integer :: muq, kq, lq
  character(4) :: senbu

  integer :: MaxNumVHS

  real(8) :: KappaEnv, pdt
  logical :: flagm

  integer :: ik, nzk
  real(8) :: stepsize
  character(3) :: cstate, skapp
  character(50) :: fnzk

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  call SRArguments(my_rank,p,tubetype,subband,kapp)

  !--------------------------
  ! START : INPUT PARAMETER
  !--------------------------

  u0 = 11.3d0
  ! the energy cost to place two electrons on a single site
  ! see eq.(14), J. Jiang et al, Phys. Rev. B 75, 035407 (2007).

  ! nanotube type
  ! s1 : Semiconductor Type I
  ! s2 : Semiconductor Type II
  ! s0 : Metal
  if(tubetype.eq.'s1') tubetypeindex = 1
  if(tubetype.eq.'s2') tubetypeindex = 2
  if(tubetype.eq.'s0') tubetypeindex = 0

  dtmin = 0.5D0
  ! dtmin : minimum diameter [nm]
  dtmax = 3.0D0
  ! dtmax : maximum diameter [nm]
  ! diameter range depends on information of the structure optimization by ETB.
  ! this information is stored in "strrel.dat".
  ! for other nanotube diameters, you run relaxte.f90 and get information.

  ! spin type
  ! s : singlet
  ! t : triplet
  stype = 's'

  ! orbital type
  ! go : A2
  ! ge : A1
  ! k1 : E
  ! k2 : E*
  otype = 'go'

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

  call itoc(int(10*kapp),skapp)
  fn3 =  'data-qex/qexkat-'//tubetype//'-e'//char(cline+AsciiZero)//&
       &char(vline+AsciiZero)//'-'//skapp//'.dat'
  fn4 =  'data-kapp/kapp-'//tubetype//'-e'//char(cline+AsciiZero)//&
       &char(vline+AsciiZero)//'-'//skapp//'.dat'
  ! output filename

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
        write(*,401) fn10
        stop
     endif

     do i=1,5
        read(21,100,iostat=io)
        if (io.ne.0) then
           write(*,402) fn10
           stop
        endif
     enddo

     do
        read(21,115,iostat=io) nn, mm, lini, tt, dt, theta, tab,&
             & phiab, dab, en, toler, iterat
        if( io.ne.0 ) exit

        if (mod(2*nn+mm,3).ne.tubetypeindex) cycle

        if( (otype=='go').or.(otype=='ge') ) then
           MaxNumVHS = numvhs(nn,mm)
           if( MaxNumVHS < subband ) cycle
        end if
        
        !call readmaruyama(nn,mm,flagm)
        !call readhata(nn,mm,flagm)
        !if( .not.flagm ) cycle
        !if( .not.((nn==6).and.(mm==6)) ) cycle

        dtinnm=dt
        dtina=dt/auc

        if (dtinnm.lt.dtmin.or.dtinnm.gt.dtmax) cycle

        ! initial data for ETB wavefunction
        call tubpar (nn,mm,length,ttini,dtini,theta,tabini,phiabini,&
             dabini,a1t,a1phi,a2t,a2phi,nc,nt,kcx,kcy,ktx,kty)
        call tbtini (tt,dt,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nab,&
             rab,hab,oab,vab)

        rank=rank+1
        if(rank.ge.p) rank=1
        if (my_rank.eq.rank) then

           kbase1=(/kcx,kcy/)
           kbase2=(/ktx,kty/)

           ! tube coordinate
           call tbcoord
           ! cutting lines around K
           call cutline(1,nmk1)
           ! cutting lines around K'
           call cutline(2,nmk2)
           ! cutting lines around G
           call qline(nmk1,nmq)

           muc = nmk1%mu(cline) - nmk1%mu(vline)

           call selectotype(otype,type)

           ! find center-of-mass momentum kc-kv
           call findkc(type,nmk1,nmk2,ao,cline,vline,nmkc)
           ! find relative motion momentum kc+kv
           call findkr(type,nmk1,nmk2,muc,kc,ao,cline,vline,nmkr)

           ! dielectric function (con, con**)
           ! con:QP con11,con12,con21:CB
           ! 11:K-->K or K'-->K', 12:K-->K', 21:K'-->K
           ! CB potential FT (vq, vq**)
           ! vq:QP vq11,vq12,vq21:CB
           ! 11:K-->K or K'-->K', 12:K-->K', 21:K'-->K
           call ftdia(muc,kc,nmkr,nmk1,nmk2,nmq,nmwf2,nmwf3,&
                &wfcoe1,wfcoe2,wfcoe3,con,con11,con12,&
                &con21,vq,vqc11,vqc12,vqc21)

           ! write dielectric function
           senbu =  num2str(nn,mm)
           fn5 = 'data-pol/polar'//'-'//senbu//'-'//skapp//'.dat'
           open(20,FILE=fn5)
           do lq = 1, nmq%nline
              muq = nmq%mu(lq)
              if( muq /= 0 ) cycle
              !write(20,*) muq, nmq%bt, nmq%up
              !write(20,*) 1.0D0/dt, 1.0D0+con(muq,kq)/kapp
              do kq = nmq%bt, nmq%up
                 write(20,*) kq*((2.0d0*pi/tt)/nt), 1.0D0+con(muq,kq)
              end do
              write(20,*)
           end do
           close(20)


           if(type.eq.'ga') then
              ! A excitons
              ! ec:conduction QP, ev:valence QP, ecs:conduction SP
              ! ecv:valence SP
              ! wc:Kd, wx:Kx
              call iniga(muc,kc,nmkr,nmq,con,con11,con12,&
                   &con21,vq,vqc11,vqc12,vqc21,nmwf2,nmwf3,&
                   &wfcoe1,wfcoe2,wfcoe3,ec,ev,ecs,evs,wc,wx,kcv,fcv)
           else if( (type=='k1').or.(type=='k2') ) then
              call inik(muc,kc,nmkr,ec,ev,wc,wx)
           end if

           nd=nmukr

           ! exciton eigenvalues and eigenvectors
           ! ex(state):velues exvec(:,state):vectors
           call exh(otype,stype,nd,nmkr,ec,ev,wc,wx,ex,exvec)

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

           !write exciton wave function coefficient
           stepsize = (2.0D0*PI/ttini)/(nt*5)
           call itoc(state,cstate)
           !fnzk = 'zk-swnt'//num2str(nn,mm)//'-state'//cstate//'.dat'
           fnzk = 'data-wf/zk-e'//char(cline+AsciiZero)//&
                char(vline+AsciiZero)//'-'//senbu//'-'//skapp//'.dat'
           nzk = 100*nn + mm
           open(nzk,FILE=fnzk)
           !write(nzk,FMT='(A1,2I5)') '#', nn, mm
           !write(nzk,FMT='(A1,F12.5)') '#', kapp
           do ik = 1, nd
              write(nzk,FMT='(3F15.8)') (kcv(ik,2,1,1)-nmk1%vhsk(cline))&
                   &*stepsize,abs(dble(exvec(ik,state)))
           end do
           close(nzk)

           call dgaussls(fnzk,lk)
           call findcon(fn5,lk,polareps)

           !dt
           mxs(1) = dtinnm
           !theta
           mxs(2) = theta
           !excitation energy
           mxs(3) = ex(state)
           !binding energy
           mxs(4) = ecvc-ex(state)
           !QP energy
           mxs(5) = ecvc
           !selfenergy
           mxs(6) = ecvc-ecvsc
           !selfenergy-bd
           mxs(7) = mxs(6)-mxs(4)
           !exciton size in k-space
           mxs(8) = lk
           !polarization function at lk
           mxs(9) = polareps
           !q vector
           qvecs(:) = nmk1%mu(subband)*kbase1(:) + nmk1%vhsk(subband)*kbase2(:)

           dest=0
           call MPI_SEND (mxs,size,MPI_COMPLEX16,dest,tag,MPI_COMM_WORLD,ierr)
           call MPI_SEND (qvecs,sizeq,MPI_REAL8,dest,tagq,MPI_COMM_WORLD,ierr)

           call MPI_SEND(kapp,sizek,MPI_REAL8,dest,tagk,MPI_COMM_WORLD,ierr)

        end if

     end do

     close (21,status='keep')

  end if


  if( my_rank.eq.0 ) then
     ! gather results from other machines and write results to output file

     open(53,file=fn3,status='unknown')
     close(53,status='delete')
     open(54,file=fn4,status='unknown')
     close(54,status='delete')

     open (53,file=fn3,status='new')
     open (54,file=fn4,status='new')

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
        read(21,115,iostat=io) nn,mm,lini,tt,dt,theta,tab,phiab,&
             &dab,en,toler,iterat
        if (io.ne.0) exit

        if (mod(2*nn+mm,3).ne.tubetypeindex) cycle

        if( (otype=='go').or.(otype=='ge') ) then
           MaxNumVHS = numvhs(nn,mm)
           if( MaxNumVHS < subband ) cycle
        end if

        !call readhata(nn,mm,flagm)
        !if( .not.flagm ) cycle
        !if( .not.((nn==6).and.(mm==6)) ) cycle
        
        dtinnm = dt
        dtina = dt / auc
        if(dtinnm.lt.dtmin.or.dtinnm.gt.dtmax) cycle

        rank = rank + 1
        if(rank.ge.p) rank = 1
        source = rank
        call MPI_RECV (mxr,size,MPI_COMPLEX16,source,tag,&
             &MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_COMPLEX16, count, ierr )
        call MPI_RECV (qvecr,sizeq,MPI_REAL8,source,tagq,&
             &MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_REAL8, countq, ierr )

        call MPI_RECV (kapp,sizek,MPI_REAL8,source,tagk,&
             &MPI_COMM_WORLD,status,ierr)
        call MPI_Get_count ( status, MPI_REAL8, countk, ierr )

        write(53,110) nn, mm, mxr(:), kapp, qvecr(:)
        write(54,111) nn, mm, mxr(1:3), mxr(8), kapp
        ! write calculated results to the file named ``fn3''
        ! see INPUT in the head of this source code.

        write(*,FMT='(A10,2I4)') "finish: ", nn, mm
     end do

     close (21,status='keep')

  end if


  call MPI_FINALIZE (ierr)


110 format(2i4,9f10.4,f7.1,2f10.4)
111 format(2i4,5f10.4,f7.1)
120 format(6f12.5)
115 format(" ",ss,i3.3," ",ss,i3.3," ",sp,d22.15," ",sp,d22.15," ",sp,d2&
       &2.15," ",sp,d22.15," ",sp,d22.15," ",sp,d22.15," ",sp,d22.15," ",s&
       &p,d22.15," ",sp,d22.15," ",ss,i6.6," ")
100 format(a223)
401 format(" file ",a10," is not found >> ")
402 format(" file ",a10," has wrong format >> ")
502 format(" (",ss,i3.3,",",ss,i3.3,") ")

close(53)
close(54)

end program main
