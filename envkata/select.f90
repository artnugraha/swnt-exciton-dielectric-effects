!>
!! \file select.f90
!<

!>
!! \brief select exciton type
!!
!! \param[in] otype orbital type
!! \param[out] exciton type
!<
subroutine selectotype(otype,type)

  implicit none

  character(len=2),intent(in) :: otype
  character(len=2),intent(out) :: type

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

end subroutine selectotype
