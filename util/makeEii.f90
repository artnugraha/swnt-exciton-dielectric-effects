MODULE dbase
! single database variable
  ! Eii(n,m,i,nkapp,flagEii)
  ! dttheta(n,m,flagdt)

  REAL(8) :: Eii(40,0:30,4,80,2)
  ! Eii(n,m,i,nkapp,flagEii)
  ! flagEii = 1 --> Eii value
  ! flagEii = 2 --> lk value

  ! for metallic SWNT:
  ! i = 1 --> M11L, i = 2 M11H

  ! for semiconducting SWNT:
  ! i = 1 --> S11, i = 2 --> S22, i = 3 --> S33, i = 4 --> S44

  REAL(8) :: dttheta(40,0:30,2)
  ! dttheta(n,m,flagdt)
  ! flagdt = 1 --> dt value
  ! flagdt = 2 --> theta value

  INTEGER :: flag(40,0:30,4,80,2)
  ! This is used for safety confirmation

  INTEGER :: n, m, ii, nkapp, ntype, j
  REAL(8) :: dt, theta, E, lk, kapp

END MODULE DBASE

PROGRAM makeEii
! make a single database for Eii

  USE dbase
  IMPLICIT NONE

!==================================================================!
! VARIABLES DECLARATION
!------------------------------------------------------------------!
! general variables:

  INTEGER :: i, dummy
 
!------------------------------------------------------------------!
! variables for reading sx-eXX.dat:

  CHARACTER(LEN=25) :: filetoread
  CHARACTER(15), PARAMETER :: root = "kappa/"

!==================================================================!
! DATA PROCESSING
!------------------------------------------------------------------!
! initialization:
  Eii = 0.d0; dttheta = 0.d0 ! gives 0 values to the dbase arrays
  flag = 0
!------------------------------------------------------------------!
! semiconducting SWNT
  ! type I
  DO i = 1, 4
     filetoread = TRIM(root)//"s1-e"//CHAR(48+i)//&
          &CHAR(48+i)//".dat"
     CALL ReadWrite(filetoread,i)
  END DO

  ! type II
  DO i = 1, 4
     filetoread = TRIM(root)//"s2-e"//CHAR(48+i)//&
          &CHAR(48+i)//".dat"
     CALL ReadWrite(filetoread,i)
  END DO

!------------------------------------------------------------------!
! metallic SWNT
  DO i = 2, 3
     filetoread = TRIM(root)//"s0-e"//CHAR(48+i)//&
          &CHAR(48+i)//".dat"
     CALL ReadWrite(filetoread,i-1)
  END DO
!------------------------------------------------------------------!
! WRITE TO OUTPUT FILE

  OPEN(UNIT=60,FILE="eii.dat",FORM="unformatted")
  WRITE(60) Eii, dttheta
  CLOSE(60)

END PROGRAM makeEii

!==================================================================!

SUBROUTINE ReadWrite(filer,itrans)
  
  USE dbase
  IMPLICIT NONE
  INTEGER :: nline ! number of lines in a file
  INTEGER, INTENT(IN) :: itrans ! transition label
  CHARACTER(LEN=25) :: filer ! filename

  OPEN(UNIT=10,FILE=TRIM(filer),STATUS="old")
  READ(10,*) ntype, ii, nline ! 1st line -> info on the source
  PRINT*, filer, nline, ntype, ii
  IF ((ii /= itrans) ) &
       STOP "Inconsistent data!"
  DO j = 1, nline
     READ(10,*), n, m, dt, theta, E, lk, kapp
     nkapp = int(10.d0 * (kapp - 1.d0 + 0.0001d0) + 1.d0)  
     ! safety confimation
     IF (n > 40) STOP "n is greater than 40"
     IF (m > 30) STOP "m is greater than 30"
     IF (flag(n,m,itrans,nkapp,1) == 0) THEN
        flag(n,m,itrans,nkapp,1) = 1
     ELSE
        STOP "Something wrong with variables assignment"
     END IF
     Eii(n,m,itrans,nkapp,1) = E
     Eii(n,m,itrans,nkapp,2) = lk
     dttheta(n,m,1) = dt
     dttheta(n,m,2) = theta
  END DO
  CLOSE(10)

END SUBROUTINE ReadWrite

!==================================================================!
