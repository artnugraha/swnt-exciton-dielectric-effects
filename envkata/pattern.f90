
program main

  implicit none

  call MakeFamilyPattern()

end program main


subroutine MakeFamilyPattern()

  implicit none

  integer,parameter :: max = 1000
  integer :: nm(2,max), n, m, fp, line, i, old
  real(8) :: mx(7,max), vec(3,max), kappa(max)

  mx(:,:) = -1.0D0
  nm(:,:) = -1
  do i = 1, max
     read(*,*,END=1000) fp, nm(:,i), mx(:,i), vec(:,i), kappa(i)
  end do

1000 line = i - 1

  old = -1
  do i = 1, line
     fp = 2*nm(1,i) + nm(2,i)
     if( old /= fp ) write(*,*) ""
     write(*,FMT=2000) fp, nm(:,i), mx(1,i), mx(2:7,i), vec(:,i), kappa(i)
     old = fp
  end do

2000 FORMAT(3I3,11F9.5)

end subroutine MakeFamilyPattern
