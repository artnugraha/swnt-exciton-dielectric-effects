
!= qline.f90
! * make cuttinge line around Gamma point
!
!=== Notice
! * module common must be complied first
!
!=== Input
! * nmk1: the cutting line information around K point from subroutine cutline
!
!=== Output
! * nmq
! ** nmq%mu(i) : cutting line number (mu)
! ** nmq%bt(i) : k_vector at starting point on the cutting line
! ** nmq%up(i) : k_vector at ending point on the cutting line
! ** nmq%nline : the number of cuttinge line (number of array)
!
!=== Contact
! author:: Jie Jiang
! e-mail:: mailto:jiang@chips.ncsu.edu
!
! modefy and comment:: Kentaro Sato
! e-mail:: mailto:kentaro@flex.phys.tohoku.ac.jp
subroutine qline(nmk1,nmq)

  use common,only:nmmesh,qmesh

  implicit none

  ! input
  type(nmmesh),intent(in) :: nmk1

  ! output
  type(qmesh),intent(out) :: nmq

  ! local
  integer :: nline, bzbt, bzup
  integer,pointer :: bzbtarr(:), bzuparr(:)
  integer :: i

  nline = nmk1%nline
  allocate(bzbtarr(nline),bzuparr(nline))

  do i = 1, nline
     bzbtarr(i) = nmk1%bt(i)
     bzuparr(i) = nmk1%up(i)
  end do

  bzbt = MINVAL(bzbtarr)
  bzup = MAXVAL(bzuparr)

  nmq%bt = -(bzup-bzbt)
  nmq%up = bzup - bzbt

  if( mod(nline,2).eq.1 ) then
     nmq%nline = nline
     do i = 1, nmq%nline
        nmq%mu(i) = -(nline-1)/2 + (i-1)
     end do
  else
     nmq%nline = nline + 1
     do i = 1, nmq%nline
        nmq%mu(i) = -nline/2 + (i-1)
     end do
  end if

  deallocate(bzbtarr,bzuparr)

end subroutine qline
