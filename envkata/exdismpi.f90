
!>
!! \file exkatampi.f90
!>

!>
!! \mainpage
!!
!! \section Input
!! -# Input Parameters
!!  - variables for exkatmpi system.
!!    See and Check Input Parameters Section in exkatampi.f90.
!! -# strrel.dat
!!  - parameter file for ETB. 
!!    you must use this, otherwise you cannot get family behaivour. 
!!    Please see Geroge-san's PhD thesis.
!! -# exkataparameter.f90
!!  - parameters for exkatmpi system.
!!    Probably you need not change parameters in exkataparameter.f90.
!!
!! \section Output
!! Output file is qex-kat-sX-eXX-XXX.dat
!! - Filename
!!  - sX is 
!!   - Metal : s0
!!   - Type I : s1
!!   - Type II : s2
!!  - eXX is Eii
!!  - XXX is exciton type
!!   - A1 : ge
!!   - A2 : go
!! - Meaing of each row of qex-kat-sX-eXX-XXX.dat
!!  -# n
!!  -# m
!!  -# diameter [nm]
!!  -# chiral angle [rad]
!!  -# exciton energy [eV] \f$ \Omega _n \f$
!!  -# binding energy [eV] \f$ E_{\rm bd} = E_{\rm QP} - \Omega _0 \f$
!!  -# quasi particle energy [eV] \f$ E({\bf k}_{\rm c}) = \epsilon({\bf k}_{\rm c}) + \Sigma({\bf k}_{\rm c}) \f$
!!  -# self energy [eV] \f$ \Sigma = \sum W \f$
!!  -# self energy - binding energy [eV] \f$ E_{\rm mb} = \Sigma - E_{\rm bd} \f$
!!
!! \attention
!! -# module common must be complied first
!! -# This program is based on MGM a1176. However some equation is not 
!!   correct. Please be careful. KS wrote comment in this source code or 
!!   in PhD thesis.
!! -# Original programs have been developped by using Intel Fortran Compiler (IFC) 7.1.
!!    Probably you need IFC 7.1...
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
!<

!>
!! \brief main program for Exciton Energy Kataura Plot
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
program main

  use common
  use exkataif,only : cutline, qline, ftdia, iniga, exh, findkc, findkr, num2str, numvhs, inik
  use exkataparameter,only : PI, AUC

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
  character(len=2) :: tubetype

  !> How many cutting lines do you take into calculation?
  !! Generary one cutting line is enough.
  !! - one : one cutting lines
  !! - all : all cutting lines
  character(len=3) :: ao

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

  !> maximau diameter [nm]
  real(8) :: dtmin
  !> minimum diameter [nm]
  real(8) :: dtmax

  !> fn3  : output file name (qex-kata-...)
  character(len=100) :: fn3
  !! fn10 : input parameter file name (strrel.dat)
  character(len=100) :: fn10
  character(len=100) :: fn20
  ! fn3  : output file name
  ! fn10 : input file name
  ! see INPUT PARAMETER in this source code

  integer :: subband

  !
  ! Output
  !
  real(8) :: ecvc, ecvsc

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
  complex(8) :: ec(2,bound),ev(2,bound), wc(2,2,bound,bound), wx(2,2,bound,bound)
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
  integer(4),parameter :: size = 10, sizeq = 2, sizek = 1
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
  character(len=3) :: cstate
  character(len=50) :: fnzk
  character(len=4) :: skapp

  integer :: tnn, tmm
  integer :: iline, iknum, kmx, liner, kcr
  real(8) :: exs, kcs, kadd

  integer,allocatable :: kmxarr(:), knum(:)
  real(8),allocatable :: expl(:,:,:)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  !call SRArguments(my_rank,p,tubetype,subband)
  call SRArguments(my_rank,p,tubetype,subband,kapp)

  !--------------------------
  ! START : INPUT PARAMETER
  !--------------------------
  
  nn = 6
  mm = 5

  u0 = 11.3d0
  ! the energy cost to place two electrons on a single site
  ! see eq.(14), J. Jiang et al, Phys. Rev. B 75, 035407 (2007).

  !kapp = 2.22d0
  ! a static dielectric constant describing the effects of electrons in core states
  ! see eq.(6) and the following paragraph, J. Jiang et al, Phys. Rev. B 75, 035407 (2007).

  !KappaEnv = 3.0D0
  ! dielectric constant around nanotube

  ! nanotube type
  ! s1 : Semiconductor Type I
  ! s2 : Semiconductor Type II
  ! s0 : Metal
  if(tubetype.eq.'s1') tubetypeindex = 1
  if(tubetype.eq.'s2') tubetypeindex = 2
  if(tubetype.eq.'s0') tubetypeindex = 0

  dtmin = 0.6D0
  ! dtmin : minimum diameter [nm]
  dtmax = 1.6D0
  ! dtmax : maximum diameter [nm]
  ! diameter range depends on information of nanotube structure optimazation by ETB.
  ! this information is stored in "strrel.dat".
  ! if you want to calculate large diameter nanotube, you run relaxte.f90 and get information.

  ! sping type
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

  !fn3 = 'qex-kat-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
  !     &//char(state-1+AsciiZero)//'.dat'

  call rtos(kapp,skapp)
  fn3 = 'qex-dis-'//tubetype//'-e'//char(cline+AsciiZero)//char(vline+AsciiZero)//'-'//otype&
       &//char(state-1+AsciiZero)//'-kappa'//skapp//'.dat'
  ! output filename

  fn10 = './strrel.dat'
  ! input data from ETB

  !------------------------
  ! END : INPUT PARAMETER
  !------------------------

  open(21,file=fn10,status='old',iostat=io)
  if (io.ne.0) then
     write(*,401)fn10
     call MPI_FINALIZE (ierr)
     stop
  endif

  do i=1,5
     read(21,100,iostat=io)
     if (io.ne.0) then
        write(*,402)fn10
        call MPI_FINALIZE (ierr)
        stop
     endif
  enddo

  do
     read(21,115,iostat=io) tnn, tmm, lini, tt, dt, theta, tab, phiab, dab, en, toler, iterat
     if( io.ne.0 ) then
        call MPI_FINALIZE (ierr)
        stop
     end if
     
     if( (tnn==nn).and.(tmm==mm) ) exit
  end do

  close(21)
  
  if (mod(2*nn+mm,3).ne.tubetypeindex) then
     call MPI_FINALIZE (ierr)
     stop
  end if
  
  if( (otype=='go').or.(otype=='ge') ) then
     MaxNumVHS = numvhs(nn,mm)
     if( MaxNumVHS < subband ) then
        call MPI_FINALIZE (ierr)
        stop
     end if
  end if
        
  dtinnm=dt
  dtina=dt/auc

  if (dtinnm.lt.dtmin.or.dtinnm.gt.dtmax) then
     call MPI_FINALIZE (ierr)
     stop
  end if
  
  ! initial data for ETB wavefunction
  call tubpar (nn,mm,length,ttini,dtini,theta,tabini,phiabini,&
       dabini,a1t,a1phi,a2t,a2phi,nc,nt,kcx,kcy,ktx,kty)
  call tbtini (tt,dt,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nab,&
       rab,hab,oab,vab)

  !call SelectKappa(tubetype,cline,dt,kapp)

  stepsize = (2.0D0*pi/tt) / dble(nt)
  kbase1=(/kcx,kcy/)
  kbase2=(/ktx,kty/)

  if( my_rank == 0 ) write(*,FMT='(A12,F15.8)') "step size =", stepsize

  ! tube coordinate
  call tbcoord
  ! cutting lines around K
  call cutline(1,nmk1)
  ! cutting lines around K'
  call cutline(2,nmk2)
  ! cutting lines around G
  call qline(nmk1,nmq)

  !muc = nmk1%mu(cline) - nmk1%mu(vline)

  call selectotype(otype,type)

  ! find center-of-mass momentum kc-kv
  call findkc(type,nmk1,nmk2,ao,cline,vline,nmkc)

  if( my_rank == 0 ) then
     write(*,FMT='(2I6)') nt, nt

     write(*,*) "around K point"
     do iline = 1, nmk1%nline
        write(*,FMT='(3I6)') nmk1%mu(iline), nmk1%bt(iline), nmk1%up(iline)
     end do

     write(*,*) "around K' point"
     do iline = 1, nmk2%nline
        write(*,FMT='(3I6)') nmk2%mu(iline), nmk2%bt(iline), nmk2%up(iline)
     end do

     write(*,*) "around Gamma point"
     do iline = 1, nmkc%nline
        write(*,FMT='(3I6)') nmkc%mu(iline), nmkc%bt(iline), nmkc%up(iline)
     end do
  end if
  
  rank = 0

  write(*,FMT='(A30,I3)') "Initialzing is finished: CPU ", my_rank

  do iline = 1, nmkc%nline
     muc = nmkc%mu(iline) ! center of mass momentum
     !muc = nmk1%mu(cline) - nmk1%mu(vline)
     kcs = 0.0D0

     !write(*,*) nmkc%nline, nmkc%bt(iline), nmkc%up(iline)

     do kc = nmkc%bt(iline), nmkc%up(iline)

        if( (kc<0).or.(abs(kc*stepsize)>4.0D0) ) cycle

        kcs = kcs + 1.0D0
        exs = 0.0D0
        rank = rank + 1
        if( rank >= p ) rank = 1

        if( my_rank == rank ) then
           ! find relative motion momentum kc+kv
           call findkr(type,nmk1,nmk2,muc,kc,ao,cline,vline,nmkr)

           ! dielectric function (con, con**)
           ! con:QP con11,con12,con21:CB
           ! 11:K-->K or K'-->K', 12:K-->K', 21:K'-->K
           ! CB potential FT (vq, vq**)
           ! vq:QP vq11,vq12,vq21:CB
           ! 11:K-->K or K'-->K', 12:K-->K', 21:K'-->K
           call ftdia(muc,kc,nmkr,nmk1,nmk2,nmq,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,con,con11,con12,con21,vq,vqc11,vqc12,vqc21)

           if(type.eq.'ga') then
              ! A excitons
              ! ec:conduction QP, ev:valence QP, ecs:conduction SP
              ! ecv:valence SP
              ! wc:Kd, wx:Kx
              call iniga(muc,kc,nmkr,nmq,con,con11,con12,con21,vq,vqc11,vqc12,vqc21,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,ec,ev,ecs,evs,wc,wx,kcv,fcv)
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

           ! write exciton wave function coefficient
           !stepsize = (2.0D0*PI/ttini) / nt
           !call itoc(state,cstate)
           !fnzk = 'zk-swnt'//num2str(nn,mm)//'-kc'//'-'//otype//'-state'//cstate//'.dat'
           !nzk = 100*nn + mm
           !open(nzk,FILE=fnzk)
           !do ik = 1, nd
           !   write(nzk,FMT='(I8,2F22.15)') kcv(ik,2,1,1), dble(exvec(ik,state)), dimag(exvec(ik,state))
           !end do
           !close(nzk)

           mxs(1) = dble(iline)
           mxs(2) = kcs
           mxs(3) = kc * stepsize

           ! binding energy
           mxs(4) = ecvc - ex(state)
           ! QP energy
           mxs(5) = ecvc
           ! self energy
           mxs(6) = ecvc - ecvsc
           ! selfenergy - bd
           mxs(7) = ex(state) - ecvsc
           ! one-particle energy
           mxs(8) = ecvsc 

           ! exciton energy
           mxs(9:10) = ex(1:2)

           write(*,FMT='(2F15.8)') kc*stepsize, ex(state)

           dest=0
           call MPI_SEND(mxs,size,MPI_REAL8,dest,tag,MPI_COMM_WORLD,ierr)
           call MPI_SEND(kapp,sizek,MPI_REAL8,dest,tagk,MPI_COMM_WORLD,ierr)

        end if
     end do
  end do

  if( my_rank == 0 ) then
     ! gather results from other machines and write results to output file

     kadd = 0.0D0

     open(53,file=fn3,status='unknown')
     close(53,status='delete')
     open (53,file=fn3,status='new')

     rank=0

     allocate( kmxarr(nmkc%nline), knum(nmkc%nline) )

     do iline = 1, nmkc%nline
        kmxarr(iline) = nmkc%up(iline) - nmkc%bt(iline) + 1
     end do
     kmx = maxval(kmxarr)

     allocate( expl(nmkc%nline,kmx,size-2) )
     
     do iline = 1, nmkc%nline
        knum(iline) = 0
        do kc = nmkc%bt(iline), nmkc%up(iline)
           if( (kc<0).or.(abs(kc*stepsize)>4.D0) ) cycle

           knum(iline) = knum(iline) + 1

           rank = rank + 1
           if(rank.ge.p) rank = 1

           source = rank

           call MPI_RECV (mxr,size,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
           call MPI_Get_count ( status, MPI_REAL8, count, ierr )
           call MPI_RECV (kapp,sizek,MPI_REAL8,source,tagk,MPI_COMM_WORLD,status,ierr)
           call MPI_Get_count ( status, MPI_REAL8, countk, ierr )

           liner = idint(mxr(1))
           kcr = idint(mxr(2))

           expl(liner,kcr,:) = mxr(3:10)
           expl(liner,kcr,1) = expl(liner,kcr,1)
        end do
     end do

     do iline = 1, nmkc%nline
        do iknum = 1, knum(iline)
           write(53,FMT='(8F15.8)') expl(iline,iknum,1:8)
        end do
        write(53,*)
     end do

     deallocate (kmxarr,knum,expl)
     close(53)
  end if

  call MPI_FINALIZE (ierr)


110 format(2i5,9f12.5,f12.5)
120 format(6f12.5)
115 format(" ",ss,i3.3," ",ss,i3.3," ",sp,d22.15," ",sp,d22.15," ",sp,d2&
       &2.15," ",sp,d22.15," ",sp,d22.15," ",sp,d22.15," ",sp,d22.15," ",s&
       &p,d22.15," ",sp,d22.15," ",ss,i6.6," ")
100 format(a223)
401 format(" file ",a10," is not found >> ")
402 format(" file ",a10," has wrong format >> ")
502 format(" (",ss,i3.3,",",ss,i3.3,") ")

end program main
