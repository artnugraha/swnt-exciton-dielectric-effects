
program main

  implicit none

  call SortSB()

end program main


subroutine SortSB()

  implicit none

  integer,parameter :: max = 3000
  integer :: nm(3,max), line, i, file
  real(8) :: mx(7,max), qvec(3,max)

  mx(:,:) = -1.0D0
  nm(:,:) = -1
  do i = 1, max
     read(*,*,END=1000) nm(:,i), mx(:,i), qvec(:,i)
  end do

1000 line = i - 1

  do i = 1, line
     if( nm(1,i) < 0 ) exit
     file = 100*nm(2,i) + nm(3,i)
     write(file,FMT=2000) nm(:,i), mx(:,i), qvec(:,i)
  end do

2000 FORMAT(3I3,10F10.5)

end subroutine SortSB
