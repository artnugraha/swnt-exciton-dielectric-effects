!>
!! \file
!<


	!>
	!! \ingroup tbdftsnt
	!! \brief tight-binding dispersion for the rolled up carbon nanotube
	!!
	!! \note
	!! subroutine tbtubecb based on subroutine tbtube was created by Jiang-san.
	!!
        !! \author Jie Jiang (jiang@chips.ncsu.edu)
        !! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
        !<
	subroutine tbtubecb (evflag,kc,kt,nc,nt,nab,rab,hab,oab,ener,ampl)
	implicit none
	integer*4 ne,nb,nd,nu,ns
	parameter (ne=2,nb=4,nd=3,nu=64,ns=2)
	integer*4 nm
	parameter (nm=nb*ns)
	integer*4 evflag,nc,nt
	real*8 kc,kt
!
	integer*4 nab(ns,ns)
	real*8 rab(nd,nu,ns,ns)
	real*8 hab(nb,nb,nu,ns,ns)
	real*8 oab(nb,nb,nu,ns,ns)
	real*8 ener(nm)
	complex*16 ampl(nm,nm)
	real*8 onsite(ne)
!	parameter (onsite=(/-0.13573d+02,-0.4882d+01/))
	parameter (onsite=(/-0.13573d+02,-0.5372d+01/))
	real*8 ef
	parameter (ef=-0.475044486465750d+01)
!	parameter (ef=-0.424778676d+01)
	real*8 pi
	parameter (pi=+0.3141592653589793d+01)
	complex*16 c0,c1
	parameter (c0=(+0.0d+00,+0.0d+00))
	parameter (c1=(+0.1d+01,+0.0d+00))
	integer*4 i,j,bc,bn,sc,sn,un
	real*8 fkt,fkphi
	complex*16 pf
	complex*16 h(nm,nm),o(nm,nm)
!=begin=lapack=data=!
	integer*4 n
	parameter (n=nm)
	character*1 jobz,ceval,cevec,uplo
	parameter (ceval='N',cevec='V',uplo='U')
	integer*4 itype,lda,ldb,lwork,lrwork,info
	parameter (itype=1,lda=n,ldb=n,lwork=(32+1)*n,lrwork=3*n-2)
	complex*16 a(lda,n)
	complex*16 b(ldb,n)
	real*8 w(n)
	complex*16 work(lwork)
	real*8 rwork(lrwork)
!=end=lapack=data=!


	fkt=kt/dfloat(nt)*2.0d0*pi
	fkphi=kc

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
	pf=cdexp(dcmplx(0.0d0,-fkt*rab(1,un,sn,sc)-fkphi*rab(2,un,sn,sc)))
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
	jobz=ceval
	else
	jobz=cevec
	endif
	do i=1,n
	do j=1,n
	a(j,i)=h(j,i)
	enddo
	enddo
	do i=1,n
	do j=1,n
	b(j,i)=o(j,i)
	enddo
	enddo
	call zhegv (itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
	if (info.ne.0) then
	write(*,100)
	stop
	endif
	do i=1,n
	ener(i)=w(i)-ef
	enddo
	if (evflag.ne.0) then
	do i=1,n
	do j=1,n
	ampl(j,i)=a(j,i)
	enddo
	enddo
	endif
!=end=lapack=code=!

	return
100	format(" >> tbdftsnt >> tbtube >> lapack >> zhegv >> info >> overflo&
     &w >> ")

	end subroutine tbtubecb



