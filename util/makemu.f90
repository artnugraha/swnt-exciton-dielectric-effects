PROGRAM makemu

  CHARACTER(LEN=*), PARAMETER :: fni = "efmrel.dat"
  CHARACTER(LEN=*), PARAMETER :: fno = "mu.dat"
  INTEGER :: n, m, i, iread
  REAL(8) :: masse, massh, masseh
  REAL(8) :: mu(40,0:30,20,3)

  ! initialization
  mu = 0.0d0
  OPEN(UNIT=10,FILE=fni,STATUS="old",IOSTAT=io)
  DO iread = 1, 5705 ! total lines in "efmrel.dat"
     READ(10,*) n, m, i, masse, massh, masseh
     mu(n,m,i,1) = masse
     mu(n,m,i,2) = massh
     mu(n,m,i,3) = masseh
  END DO
  CLOSE(10)

  OPEN(UNIT=20,FILE=fno,FORM="unformatted")
  WRITE(20) mu
  CLOSE(20)

  ! recovery
  mu = 0.0d0
  OPEN(UNIT=30,FILE=fno,STATUS="old",FORM="unformatted")
  READ(30) mu
  CLOSE(30)

  ! check
  n = 8; m = 8; i = 2
  WRITE(*,"(3i4,3d15.5)") n, m, i, mu(n,m,i,1), mu(n,m,i,2), mu(n,m,i,3)

END PROGRAM makemu
