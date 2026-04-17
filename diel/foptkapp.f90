MODULE vardat

  ! Parameters
  INTEGER, PARAMETER :: nkapp = 80, max = 300
  REAL(8), PARAMETER :: hbar = 6.58211899d-16 ! (in eV s)
  REAL(8), PARAMETER :: me = 0.510998910d6 ! (in eV / c^2)
  REAL(8), PARAMETER :: kB = 0.861734315d-5 ! (in eV / K)
  REAL(8), PARAMETER :: c0 = 2.99792458d8 ! (in m/s)
  REAL(8), PARAMETER :: Lt = 200d-9 ! (in nm)
  REAL(8), PARAMETER :: Tp = 300 ! (room temperature in kelvin)
  !REAL(8), PARAMETER :: Mp = 1.0d-4 ! (total multiplier)
  REAL(8) :: pi, cF

  REAL(8) :: Eii(40,0:30,4,nkapp,2), dttheta(40,0:30,2)
  REAL(8) :: mu(40,0:30,20,3)
  CONTAINS

    SUBROUTINE init()
      OPEN(UNIT=60,FILE="eii.dat",STATUS="old",FORM="unformatted")
      READ(60) Eii, dttheta
      CLOSE(60)
  
      OPEN(UNIT=65,FILE="mu.dat",STATUS="old",FORM="unformatted")
      READ(65) mu
      CLOSE(65)
    END SUBROUTINE init

END MODULE vardat

PROGRAM foptkapp
! find optimized kappa from experimental Eii value

  USE vardat
  IMPLICIT NONE

!==================================================================!
! VARIABLES DECLARATION
!------------------------------------------------------------------!
  ! variables for reading data files
  INTEGER :: n(max), m(max), p(max), i(max), ndata, nselect
  REAL(8) :: dt(max), Ecal(max), Eexp(max), lk(max), xvar
  CHARACTER(LEN=8) :: input
  CHARACTER(LEN=50) :: comment

  ! variables for data processing
  INTEGER :: ikapp, j, nopt, nopt1, nopt2, nsign(max), moveout
  REAL(8) :: kapp1, kapp2, kapp(nkapp), diff(nkapp), E1, E2
  REAL(8) :: optkapp(max), deltaE, mass
  pi = 4.0d0 * DATAN(1.0d0)
  cF = ( Lt/(hbar*c0) ) * DSQRT((kB * Tp)**3) * DSQRT(me/(8*pi))  

!==================================================================!
! READ FILES
!------------------------------------------------------------------!
! UNFORMATTED DAT (eii.dat and mu.dat)
  call init()

!------------------------------------------------------------------!

! EXPERIMENTAL DATA (hata.dat)
! format of the data is: n, m, Eii(exp), p, i
  
  PRINT "(a)",&
       "Enter a source (*.dat): "
  READ *, input
  input = trim(input)//".dat"
  OPEN(UNIT=70,FILE=input,STATUS="old")
  READ(70,"(a)") comment ! skip 1st line
  READ(70,*), ndata ! read the number of data in a file
  DO j = 1, ndata
     READ(70,*) n(j), m(j), Eexp(j), p(j), i(j)
     IF (j == 1) THEN
        mass = 1.0d0 / (1.d0/mu(n(j),m(j),i(j),1)&
             +1.d0/mu(n(j),m(j),i(j),2))
        deltaE = cF * DSQRT(mass)
        Eexp(j) = Eexp(j) + deltaE
     END IF
  END DO
  CLOSE(70)

!==================================================================!
! DATA PROCESSING
!------------------------------------------------------------------!
  
  ! change kappa integer (location in dbase) to kappa real
  DO ikapp = 1, nkapp
     kapp(ikapp) = 1.0d0 + 0.1d0 * DBLE(ikapp - 1)
  END DO
  
  nsign = 1
  moveout = 0
  dataloop: DO j = 1, ndata 
     DO ikapp = 1, nkapp
        diff(ikapp) = abs(Eexp(j) - Eii(n(j),m(j),i(j),ikapp,1))
     END DO
     nopt = MINLOC(diff,dim=1)
     diff(nopt) = Eexp(j) - Eii(n(j),m(j),i(j),nopt,1)
     IF (diff(nopt) < 0) THEN
        nopt1 = nopt
        nopt2 = nopt1 + 1
        E1 = Eii(n(j),m(j),i(j),nopt1,1)
        E2 = Eii(n(j),m(j),i(j),nopt2,1)
        kapp1 = kapp(nopt1)
        kapp2 = kapp(nopt2)
        optkapp(j) = kapp1 + (kapp2 - kapp1) * &
             ( (Eexp(j)-E1) / (E2-E1) )
     ELSE IF (diff(nopt) > 0) THEN
        nopt1 = nopt - 1
        IF (nopt1 < 1) THEN
           optkapp(j) = 1.0d0
           moveout = moveout + 1
           nsign(j) = 0
        ELSE
           nopt2 = nopt
           E1 = Eii(n(j),m(j),i(j),nopt1,1)
           E2 = Eii(n(j),m(j),i(j),nopt2,1)
           kapp1 = kapp(nopt1)
           kapp2 = kapp(nopt2)
           optkapp(j) = kapp1 + (kapp2 - kapp1) * &
                ( (Eexp(j)-E1) / (E2-E1) )           
        END IF
     ElSE
        optkapp(j) = kapp(nopt)   
     END IF
  END DO dataloop

  OPEN(71,file="cal-"//TRIM(input))
  WRITE(71,"(a)") comment
  WRITE(71,"(i5)") ndata - moveout  
  loopwrite : DO j = 1, ndata
     IF (nsign(j) == 1) THEN
        dt(j)   = dttheta(n(j),m(j),1)
        ikapp   = int(10.d0*(optkapp(j)-1.d0+0.0001d0)+1.d0)
        Ecal(j) = Eii(n(j),m(j),i(j),ikapp,1)
        lk(j)   = Eii(n(j),m(j),i(j),ikapp,2)
        xvar    = p(j)**0.8 *  (1/dt(j))**1.6 * (1/lk(j))**0.4
        WRITE(71,100) n(j), m(j), i(j), p(j), dt(j), lk(j),&
             xvar, optkapp(j), Ecal(j), Eexp(j)
        WRITE(*,200) j, n(j), m(j), p(j), i(j), Eexp(j), Ecal(j), optkapp(j)
        CYCLE loopwrite
     END IF
  END DO loopwrite

  CLOSE(71)

100 FORMAT (4i3,6f8.4)
200 FORMAT (5i3,3f10.4)
 
END PROGRAM foptkapp
