
subroutine readhata(nn,mm,flag)

  implicit none

  integer :: nn, mm
  logical :: flag
  integer :: n, m, i

  flag = .FALSE.
  open(235,FILE='hata.txt')

  do i = 1, 175
     read(235,*) n, m
     if( (nn==n).and.(mm==m) ) flag = .TRUE.
  end do
  
  close(235)

end subroutine readhata


subroutine readmaruyama(nn,mm,flag)

  implicit none

  integer :: nn, mm
  logical :: flag
  integer :: n, m, i

  flag = .FALSE.
  open(235,FILE='maruyama.txt')

  do i = 1, 82
     read(235,*) n, m
     if( (nn==n).and.(mm==m) ) flag = .TRUE.
  end do
  
  close(235)

end subroutine readmaruyama
