
program main

  implicit none

  call MakeFamily()

end program main


subroutine MakeFamily()

  implicit none

  integer,parameter :: max = 1000
  integer :: nm(2,max), n, m, fp, line, i, old
  real(8) :: mx(7,max), qvec(2,max), kappa(max)

  mx(:,:) = -1.0D0
  nm(:,:) = -1
  do i = 1, max
     read(*,*,END=1000) nm(:,i), mx(:,i), qvec(:,i), kappa(i)
  end do

1000 line = i - 1

  do i = 1, line
     if( nm(1,i) < 0 ) exit
     fp = 2*nm(1,i) + nm(2,i)
     write(*,FMT=2000) fp, nm(:,i), mx(1,i), mx(2:7,i), qvec(:,i), kappa(i), sqrt(dot_product(qvec(:,i),qvec(:,i)))
  end do

2000 FORMAT(3I3,11F9.5)

end subroutine MakeFamily
