
subroutine mukrg(nline1,mu1,bt1,up1,nline2,mu2,bt2,up2,nline,mu,bt,up)

  use common,only : mu2d, mumx
  ! STB: mu2d -> murmx

  implicit none

  ! input
  integer,intent(in) :: nline1, mu1(mu2d), bt1(mu2d), up1(mu2d)
  integer,intent(in) :: nline2, mu2(mu2d), bt2(mu2d), up2(mu2d)

  ! output
  integer,intent(out) :: nline, mu(mumx), bt(mumx), up(mumx)

  ! local
  integer,allocatable :: mucom(:)
  integer :: nmu, ncom, mut, mut1, mut2, arr(4), min, max, btt, upt
  integer :: i, j, i1

  nmu = 0
  ncom = 0

  allocate(mucom(nline1*nline2))

  do i = 1, nline1
     do j = 1, nline2
        mut = mu1(i) - mu2(j)
        ncom = ncom + 1
        mucom(ncom) = mut
        if( ncom.eq.1 ) then
           nmu = nmu + 1
           mu(nmu) = mut
        else
           do i1 = 1, ncom - 1
              if( mut.eq.mucom(i1) ) goto 10
           end do
           nmu = nmu + 1
           mu(nmu) = mut
        endif
        if (nmu.gt.mumx) then
           write (*,100)
           stop
        end if
10   end do
  end do

  nline=nmu

  do i1=1,nmu
     mut1=mu(i1)
     min=10000
     max=-10000
     do i=1,nline1
        do j=1,nline2
           mut2=mu1(i)-mu2(j)
           if (mut2.eq.mut1) then
              arr(1)=bt1(i)-bt2(j)
              arr(2)=bt1(i)-up2(j)
              arr(3)=up1(i)-bt2(j)
              arr(4)=up1(i)-up2(j)
              btt=MINVAL(arr)
              upt=MAXVAL(arr)
              if (btt.le.min) min=btt
              if(upt.gt.max) max=upt
           endif
        enddo
     enddo
     bt(i1)=min
     up(i1)=max
  enddo

100 format(" >> mukrg >> mumx too small >> ")

  deallocate (mucom)

end subroutine mukrg
