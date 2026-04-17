PROGRAM delta

  IMPLICIT NONE

  REAL(8), PARAMETER :: kT = 26 ! meV
  INTEGER, PARAMETER :: mesh  = 1000
  INTEGER, PARAMETER :: iter  = 100

  REAL(8) :: sumup, sumdown, x, fx
  INTEGER :: n, i
 
  OPEN(10,FILE="checkdel.dat",STATUS="replace")
  x = 0.0d0
  DO i = 1, mesh
     sumup = 0.0d0
     sumdown = 0.0d0
     DO n = 1, iter
        sumup = sumup + x * DBLE(n**2) * EXP(-x*DBLE(n**2))
        sumdown = sumdown + EXP(-x*DBLE(n**2))
     END DO
     fx = kT * (sumup / sumdown)
     WRITE(*,"(2f10.4)") x, fx
     WRITE(10,"(2f10.4,i4)") x, fx, i
     x = x + (1.d0 / DBLE(mesh))
  END DO
  CLOSE(10)
  
END PROGRAM delta
