program main

  implicit none

  real(8) :: sum, temp
  integer :: i

  sum = 0.0D0
  do i = 1, 159201
     read(*,*) temp
     sum = sum + temp
     write(*,*) sum
  end do
  write(*,*) sum

end program main
