module mdatEii

implicit none

!==================================================================!
! VARIABLES DECLARATION
!------------------------------------------------------------------!
! single database variable
! Eii(n,m,i,nkapp,flagEii)
! dttheta(n,m,flagdt)
  
  integer :: i, j, nkapp
  
  real(8) :: Eii(40,0:30,4,80,2)
  ! Eii(n,m,i,nkapp,flagEii)
  ! flagEii = 1 --> Eii value
  ! flagEii = 2 --> lk value

  ! for metallic SWNT:
  ! i = 1 --> M11L, i = 2 M11H

  ! for semiconducting SWNT:
  ! i = 1 --> S11, i = 2 --> S22, i = 3 --> S33, i = 4 --> S44

  real(8) :: dttheta(40,0:30,2)
  ! dttheta(n,m,flagdt)
  ! flagdt = 1 --> dt value
  ! flagdt = 2 --> theta value
  
!------------------------------------------------------------------!
! variables for reading sx-eXX.dat:

  integer :: nliner ! read number of lines

  integer :: n, m
  real(8) :: dt, theta, E, lk, kapp
  character(25) :: filetoread
  character(15), parameter :: root = "genkap/"

contains

  subroutine Eiidata
    ! make a single database for Eii
    implicit none

!==================================================================!
! DATA PROCESSING
!------------------------------------------------------------------!
! initialization:
    Eii = 0.d0; dttheta = 0.d0 ! gives 0 values to the dbase arrays

!------------------------------------------------------------------!
! semiconducting SWNT
    ! type I
    do i = 1, 4
       filetoread = trim(root)//"s1-e"//char(48+i)//&
            &char(48+i)//".dat"
       open(10,file=trim(filetoread),status="old")
       call numline(filetoread,nliner)
       print*, nliner        
       do j = 1, nliner
          read(10,*), n, m, dt, theta, E, lk, kapp
          if (n > 40) stop "n is greater than 40"
          if (m > 30) stop "m is greater than 30"
          nkapp = int(10 * ( kapp - 1 ) + 1)
          Eii(n,m,i,nkapp,1) = E
          Eii(n,m,i,nkapp,2) = lk
          dttheta(n,m,1) = dt
          dttheta(n,m,2) = theta
          write(*,"(3i6, 4f8.4)"), j, n, m, Eii(n,m,i,nkapp,1),&
               dttheta(n,m,1), Eii(n,m,i,nkapp,2), dttheta(n,m,2)
       end do
       close(10)
    end do
  
    ! type II
    do i = 1, 4
       filetoread = trim(root)//"s2-e"//char(48+i)//&
            &char(48+i)//".dat"
       open(10,file=trim(filetoread),status="old")
       call numline(filetoread,nliner)
       print*, nliner        
       do j = 1, nliner
          read(10,*), n, m, dt, theta, E, lk, kapp
          if (n > 40) stop "n is greater than 40"
          if (m > 30) stop "m is greater than 30"
          nkapp = int(10 * ( kapp - 1 ) + 1)
          Eii(n,m,i,nkapp,1) = E
          Eii(n,m,i,nkapp,2) = lk
          dttheta(n,m,1) = dt
          dttheta(n,m,2) = theta
          write(*,"(3i6, 4f8.4)"), j, n, m, Eii(n,m,i,nkapp,1),&
               dttheta(n,m,1), Eii(n,m,i,nkapp,2), dttheta(n,m,2)
       end do
       close(10)
    end do

!------------------------------------------------------------------!
! metallic SWNT
    do i = 1, 2
       filetoread = trim(root)//"s1-e"//char(48+i+1)//&
            &char(48+i+1)//".dat"
       open(10,file=trim(filetoread),status="old")
       call numline(filetoread,nliner)
       print*, nliner        
       do j = 1, nliner
          read(10,*), n, m, dt, theta, E, lk, kapp
          if (n > 40) stop "n is greater than 40"
          if (m > 30) stop "m is greater than 30"
          nkapp = int(10 * ( kapp - 1 ) + 1)
          Eii(n,m,i,nkapp,1) = E
          Eii(n,m,i,nkapp,2) = lk
          dttheta(n,m,1) = dt
          dttheta(n,m,2) = theta
          write(*,"(3i6, 4f8.4)"), j, n, m, Eii(n,m,i,nkapp,1),&
               dttheta(n,m,1), Eii(n,m,i,nkapp,2), dttheta(n,m,2)
       end do
       close(10)
    end do
!------------------------------------------------------------------!
! WRITE TO OUTPUT FILE

    open(60,file="eii.dat",form="unformatted")
    write(60), Eii, dttheta
    close(60)
  end subroutine Eiidata
end module mdatEii

!==================================================================!
! COUNTING NUMBER OF LINES IN A FILE
!------------------------------------------------------------------!

subroutine numline(file,nline)
  
  character(40), intent(in) :: file
  integer, intent(out) :: nline
  integer :: dummy
  open(40,file=trim(file),status="old")
  nline = 0
  do
     read(40,*,end = 400), dummy
     nline = nline + 1
  end do
400 rewind(40)
  close(40)
end subroutine numline
