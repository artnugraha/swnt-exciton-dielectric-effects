
subroutine rtos(r,s)

  implicit none

  real(8),intent(in) :: r
  character(len=4),intent(out) :: s

  integer :: i, a
  character(len=1) :: n1, n2

  a = int(r)
  n2 = char(48+a/10)
  n1 = char(48+mod(a,10))

  s = n2//n1

end subroutine rtos
