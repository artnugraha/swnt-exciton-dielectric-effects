!>
!! \file wflinecbi.f90
!! \brief calculate bounds of wave vector around K or K' point for calculation
!<

!>
!! \brief bounds of wave vector around K or K' point for calculation
!!
!! \param[in] nmk1 wave vector and energy band information around K point
!! \param[in] nmk2 wave vector and energy band information around K or K' point
!! \param[out] nmwf1 bounds of wave vector around K point for calculation
!! \param[out] nmwf2 bounds of wave vector around K or K' point for calculation
!!
!! \note
!! Since subroutine wflinecbi has the same processing,
!! subroutine wflinecbi can be divided into two parts.
!! Please try.
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
subroutine wflinecbi(nmk1,nmk2,nmwf1,nmwf2)

  use common,only : nmmesh, wfmesh

  implicit none

  ! Input
  type(nmmesh),intent(in) :: nmk1, nmk2

  ! Output
  type(wfmesh),intent(out) :: nmwf1, nmwf2

  ! Local
  integer :: nline

  nline = nmk1%nline

  nmwf1%mubt = MINVAL( nmk1%mu(1:nline) )
  nmwf1%muup = MAXVAL( nmk1%mu(1:nline) )
  nmwf1%nline = nmwf1%muup - nmwf1%mubt + 1

  nmwf1%kbt = MINVAL( nmk1%bt(1:nline) )
  nmwf1%kup = MAXVAL( nmk1%up(1:nline) )
  nmwf1%point = nmwf1%kup - nmwf1%kbt + 1



  nline = nmk2%nline

  nmwf2%mubt = MINVAL( nmk2%mu(1:nline) )
  nmwf2%muup = MAXVAL( nmk2%mu(1:nline) )
  nmwf2%nline = nmwf2%muup - nmwf2%mubt + 1

  nmwf2%kbt = MINVAL( nmk2%bt(1:nline) )
  nmwf2%kup = MAXVAL( nmk2%up(1:nline) )
  nmwf2%point = nmwf2%kup - nmwf2%kbt + 1

end subroutine wflinecbi
