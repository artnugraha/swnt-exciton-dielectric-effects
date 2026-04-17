
subroutine itoc(i,c)
  ! change data type from integer to character

  implicit none

  integer,intent(in) :: i

  integer,parameter :: MAX = 3, AsciiZero = 48
  character(len=MAX),intent(out) :: c

  integer :: j, p, q

  p = i
  do j = MAX, 1, -1
     q = mod(p,10)
     c(j:j) = achar(AsciiZero+q)
     p = p / 10
  end do

end subroutine itoc
