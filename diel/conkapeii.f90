INCLUDE "envrout.f90"

PROGRAM conkapeii

  INTEGER :: n, m, p, i, nline, j
  REAL(8) :: kappa, Eii
  CHARACTER(LEN=50) :: input, output, comment

  PRINT "(a)", "Enter an input (*.dat): "
  READ *, input
  input = trim(input)//".dat"
  OPEN(UNIT=10,FILE=input,STATUS="old")
  READ(10,"(a)") comment ! skip 1st line
  READ(10,*) nline ! read the number of data in a file
  
  PRINT "(a)", "Enter an output name (*.dat): "
  READ *, output
  output = trim(output)//".dat"
  OPEN(UNIT=20,FILE=output)
  WRITE(20,"(a)") comment
  WRITE(20,"(i5)") nline
  DO j = 1, nline
     READ(10,*) n, m, kappa, p, i
     CALL findeii(n,m,i,kappa,Eii)
     WRITE(20,"(2i3,f8.3,2i3)") n, m, Eii, p, i
  END DO
  CLOSE(10)
  CLOSE(20)
     
END PROGRAM conkapeii
