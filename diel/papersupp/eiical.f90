
INCLUDE "envrout.f90"

PROGRAM eiical
  
  USE eiimudt
  IMPLICIT NONE
  INTEGER :: n, m, i, nline, j
  REAL(8) :: fkapp, Ecal
  CHARACTER(LEN=50) :: input, output

  CALL init()

  PRINT "(a)", "Enter an input (*.dat): "
  READ *, input
  input = TRIM(input)//".dat"

  OPEN(UNIT=10,FILE=input,STATUS="old")
  READ(10,*) nline ! read the number of data in a file  

  output = "eii-"//input
  OPEN(UNIT=20,FILE=output)
  WRITE(20,"(a3,i5)") "# ", nline
  DO j = 1, nline
     READ(10,*) n, m, i, fkapp
     CALL findeii(n,m,i,fkapp,Ecal)
     WRITE(20,"(3i3,2f8.3)") n, m, i, fkapp, Ecal
  END DO
  
  CLOSE(10)
  CLOSE(20)
     
END PROGRAM eiical
