
!= findkr.f90
! * make relative motion
!
!=== Notice
! * module common must be complied first
!
!=== Input
! * type: exciton type
!   * ga: A1 and A2 exciton
!   * k1: E exciton
!   * k2: E* exciton
! * nmk1 and nmk2: cutting lines around K and K' points from subroutie cutline
! * muc: 
! * kc:
! * ao: number of calculating cutting lines
!   * all: all cutting lines
!   * one: one cutting line
!   * two: two cutting lines
! * cline and vline: number of conductance and valence bands (E_{cline}{vline})
!
!=== Output
! * nmkr
! ** nmkr%mu(i) : the relative motion of mu (nmkr%muc(i)+nmkr%muv(i))
! ** nmkr%muc(i) : mu at the conduction band
! ** nmkr%muv(i) : mu at the valence band
! ** nmkr%kr(i) : the relative motion of k (nmkr%kc(:,i)+nmkr%kv(:,i))
! ** nmkr%kc(:,i) : k at the conduction band
! ** nmkr%kv(:,i) : k at the valence band
! ** nmkr%nkr(i) : the number of 1D k point
! ** nmkr%nline : number of array
!
!=== Contact
! author:: Jie Jiang
! e-mail:: mailto:jiang@chips.ncsu.edu
!
! modefy and comment:: Kentaro Sato
! e-mail:: mailto:kentaro@flex.phys.tohoku.ac.jp

!*************************************************************************
!     find mu and k region for relative motion :
!     nline  : number of cutting lines
!     mu    : cutling line index
!     nk    : number of k point on a cutting line
!     k     : k vector
!     kc    : kc vector
!     kv    : kv vector
!*************************************************************************
subroutine findkr(type,nmk1,nmk2,muc,kc,ao,cline,vline,nmkr)

  use common,only : nmmesh, nmrel
  use common,only:mu2d,krmx,murmx
  use common,only:nmukr

  implicit none

  ! input
  integer,intent(in) :: muc, kc, cline, vline
  character(len=2),intent(in) :: type
  character(len=3),intent(in) :: ao
  type(nmmesh),intent(in) :: nmk1, nmk2

  ! output
  type(nmrel),intent(out) :: nmkr

  ! local
  integer nline1,nline2
  integer,dimension(:),allocatable::mu1,bt1,up1,mu2,bt2,up2
  integer nline,nkr,k1i,k1f,k2i,k2f,muctem,kctem,k1,k2
  integer,dimension(krmx)::kr,kcon,kval
  integer i1,i2,i,j

  if (ao.eq.'all') then

     if (type.eq.'ga') then
        nline1=nmk1%nline
        nline2=nmk2%nline
        allocate(mu1(nline1),bt1(nline1),up1(nline1),mu2(nline2),&
             bt2(nline2),up2(nline2))
        do i=1,nline1
           mu1(i)=nmk1%mu(i)
           bt1(i)=nmk1%bt(i)
           up1(i)=nmk1%up(i)
           mu2(i)=mu1(i)
           bt2(i)=bt1(i)
           up2(i)=up1(i)
        enddo
     elseif (type.eq.'k1') then
        nline1=nmk2%nline
        nline2=nmk1%nline
        allocate(mu1(nline1),bt1(nline1),up1(nline1),mu2(nline2),&
             bt2(nline2),up2(nline2))
        do i=1,nline1
           mu1(i)=nmk2%mu(i)
           bt1(i)=nmk2%bt(i)
           up1(i)=nmk2%up(i)
        enddo
        do i=1,nline2
           mu2(i)=nmk1%mu(i)
           bt2(i)=nmk1%bt(i)
           up2(i)=nmk1%up(i)
        enddo
     elseif (type.eq.'k2') then
        nline1=nmk1%nline
        nline2=nmk2%nline
        allocate(mu1(nline1),bt1(nline1),up1(nline1),mu2(nline2),&
             bt2(nline2),up2(nline2))
        do i=1,nline1
           mu1(i)=nmk1%mu(i)
           bt1(i)=nmk1%bt(i)
           up1(i)=nmk1%up(i)
        enddo
        do i=1,nline2
           mu2(i)=nmk2%mu(i)
           bt2(i)=nmk2%bt(i)
           up2(i)=nmk2%up(i)
        enddo
     endif

  endif

  select case(ao)
  case('one')
     nline1 = nmk1%nline
     nline2 = nmk2%nline
     allocate(mu1(nline1),bt1(nline1),up1(nline1),mu2(nline2),bt2(nline2),up2(nline2))
     nline1 = 1
     nline2 = 1
     if(type.eq.'ga') then
        do i=1,nline1
           mu1(i) = nmk1%mu(cline)
           bt1(i) = nmk1%bt(cline)
           up1(i) = nmk1%up(cline)
           mu2(i) = nmk1%mu(vline)
           bt2(i) = nmk1%bt(vline)
           up2(i) = nmk1%up(vline)
        enddo
     elseif (type.eq.'k1') then
        do i=1,nline1
           mu1(i)=nmk2%mu(cline)
           bt1(i)=nmk2%bt(cline)
           up1(i)=nmk2%up(cline)
        enddo
        do i=1,nline2
           mu2(i)=nmk1%mu(vline)
           bt2(i)=nmk1%bt(vline)
           up2(i)=nmk1%up(vline)
        enddo
     elseif (type.eq.'k2') then
        do i=1,nline1
           mu1(i)=nmk1%mu(cline)
           bt1(i)=nmk1%bt(cline)
           up1(i)=nmk1%up(cline)
        enddo
        do i=1,nline2
           mu2(i)=nmk2%mu(vline)
           bt2(i)=nmk2%bt(vline)
           up2(i)=nmk2%up(vline)
        enddo
     endif
  end select

  !muc=abs(muc)
  nline=0
  do i1=1,nline1
     k1i=bt1(i1)
     k1f=up1(i1)
     do i2=1,nline2
        !muctem=mu1(i1)-mu2(i2)
        !
        !muctem=abs(muctem)
        !
        !
        !if (muctem.ne.muc) cycle
        !
        k2i=bt2(i2)
        k2f=up2(i2)
        nkr=0
        do k1=k1i,k1f
           do k2=k2i,k2f
              kctem=k1-k2
              if (kctem.ne.kc) cycle
              nkr=nkr+1
              if (nkr.gt.krmx) then
                 write (*,100)
                 stop
              endif
              kr(nkr)=k1+k2
              kcon(nkr)=k1
              kval(nkr)=k2
           enddo
        enddo
        if (nkr.gt.0) then
           nline=nline+1
           if (nline.gt.murmx) then
              write (*,200)
              stop
           endif
           nmkr%mu(nline)=mu1(i1)+mu2(i2)
           nmkr%muc(nline)=mu1(i1)
           nmkr%muv(nline)=mu2(i2)
           nmkr%nk(nline)=nkr
           do i=1,nkr
              nmkr%k(nline,i)=kr(i)
              nmkr%kc(nline,i)=kcon(i)
              nmkr%kv(nline,i)=kval(i)
           enddo
        endif
     enddo
  enddo

  nmkr%nline=nline

  nmukr=0
  do i=1,nmkr%nline
     do j=1,nmkr%nk(i)
        nmukr=nmukr+1
     enddo
  enddo

  if (ao.eq.'one') then
     if (nmkr%nline.gt.1) then
        write (*,300)
        stop
     endif
  endif

100 format(" >> findkr >> nkrmx too small >> ")
200 format(" >> findkr >> nmurmx too small >> ")
300 format(" >> findkr >> nmkr%nline is not 1 >> ")

  return
end subroutine findkr
