! mueiidt <-- database for effective mass, exciton energy,
!             exciton size, swnt diameter and chiral angle

!======================================================================!
MODULE eiimudt
!======================================================================!
! database module
! see ./util directory for the raw data

  IMPLICIT NONE
  
  ! variables for reading data files
  REAL(8) :: Eii(40,0:30,4,80,2), dttheta(40,0:30,2), mu(40,0:30,20,3)
  INTEGER, PARAMETER :: nkapp = 80
  
  ! physical constants
  REAL(8), PARAMETER :: kB   = 8.61734315D-5 ! Boltzmann const. in eV/K
  REAL(8), PARAMETER :: Tr   = 300.0D0 ! room temperature in K
  REAL(8), PARAMETER :: hbar = 6.5821189916D-16 ! Planck const. in eV.s
  REAL(8), PARAMETER :: Pi   = 3.141592653589793D0
  REAL(8), PARAMETER :: c0   = 2.99792458D8 ! in m/s
  REAL(8), PARAMETER :: m0   = 0.510998910d6 ! in eV / c0^2

CONTAINS
  
  SUBROUTINE init()
    OPEN(80,FILE="eii.dat",STATUS="old",FORM="unformatted")
    READ(80) Eii, dttheta
    CLOSE(80)

    OPEN(90,FILE="mu.dat",STATUS="old",FORM="unformatted")
    READ(90) mu
    CLOSE(90)
  END SUBROUTINE init

END MODULE eiimudt
!======================================================================!
! subroutine for environmental effects calculation
!======================================================================!
SUBROUTINE findeii(nn,mm,ii,kapp,Ei)
!======================================================================!
! find Eii from kappa
  
  USE eiimudt
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nn, mm, ii
  REAL(8), INTENT(IN) :: kapp
  REAL(8), INTENT(OUT) :: Ei
  INTEGER :: ikapp

  ikapp = NINT(10.0D0 * ((kapp - 1.0D0) + 1.0D-6) + 1.D0)
  CALL INIT()
  Ei = Eii(nn,mm,ii,ikapp,1)

END SUBROUTINE findeii
!======================================================================!
SUBROUTINE findkapp(nn,mm,ii,Ei,kapp,nopt)
!======================================================================!
! find kappa from energy (rough method)

  USE eiimudt
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nn, mm, ii
  INTEGER, INTENT(OUT) :: nopt
  REAL(8), INTENT(IN)  :: Ei
  REAL(8), INTENT(OUT) :: kapp

  INTEGER :: ikapp
  REAL(8) :: diff(nkapp), kapwork(nkapp)

  CALL INIT()
  DO ikapp = 1, nkapp
     kapwork(ikapp) = 1.0D0 + 0.1D0 * DBLE(ikapp - 1)
     diff(ikapp)  = ABS(Ei - Eii(nn,mm,ii,ikapp,1))
  END DO
  nopt = MINLOC(diff,dim=1)
  kapp = kapwork(nopt)

END SUBROUTINE findkapp
!======================================================================!
SUBROUTINE foptkapp(nn,mm,ii,Ei,kapp,nopt,nout)
!======================================================================!
! find optimized kappa from energy (more accurate)

  USE eiimudt
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nn, mm, ii
  INTEGER, INTENT(OUT) :: nopt, nout
  REAL(8), INTENT(IN)  :: Ei
  REAL(8), INTENT(OUT) :: kapp

  INTEGER :: ikapp, nopt1, nopt2, isign
  REAL(8) :: diff(nkapp), kapwork(nkapp), checkdiff
  REAL(8) :: kapp1, kapp2, E1, E2

  CALL INIT()
  nout  = 0 ! number of data which will be neglected
  isign = 1
  DO ikapp = 1, nkapp
     kapwork(ikapp) = 1.0D0 + 0.1D0 * DBLE(ikapp - 1)
     diff(ikapp)  = ABS(Ei - Eii(nn,mm,ii,ikapp,1))
  END DO
  nopt = MINLOC(diff,dim=1)
  checkdiff = Ei - Eii(nn,mm,ii,nopt,1)
  IF (checkdiff < 0.0D0) THEN
     nopt1 = nopt
     nopt2 = nopt1 + 1
     E1 = Eii(nn,mm,ii,nopt1,1)
     E2 = Eii(nn,mm,ii,nopt2,1)
     kapp1 = kapwork(nopt1)
     kapp2 = kapwork(nopt2)
     kapp  = kapp1 + (kapp2 - kapp1) * &
          ( (Ei - E1) / (E2 - E1) )
  ELSE IF (checkdiff > 0.0D0) THEN
     nopt1 = nopt - 1
     IF (nopt1 == 0) THEN
        kapp  = 1.0D0
        nout  = nout + 1
        isign = 0
     ELSE
        nopt2 = nopt
        E1 = Eii(nn,mm,ii,nopt1,1)
        E2 = Eii(nn,mm,ii,nopt2,1)
        kapp1 = kapwork(nopt1)
        kapp2 = kapwork(nopt2)
        kapp  = kapp1 + (kapp2 - kapp1) * &
             ( (Ei - E1) / (E2 - E1) )
     END IF
  ELSE
     kapp = kapwork(nopt)
  END IF

END SUBROUTINE foptkapp
!======================================================================!
SUBROUTINE delta(nn,mm,x,exmass,delE)
!======================================================================!
! compute the corrections for E11S in PL spectroscopy

  USE eiimudt
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nn, mm
  INTEGER, PARAMETER   :: iter = 1000
  REAL(8), PARAMETER   :: Ln   = 2.0D-8 ! length
  REAL(8), INTENT(OUT) :: delE
  REAL(8) :: exmass, E0, x, y, sumup, sumdown
  INTEGER :: j

  CALL init()

  E0 = (hbar**2 * Pi**2 * c0**2) / (2 * m0 * Ln**2)
  y  = E0 / (kB * Tr)
  exmass = mu(nn,mm,1,3)
  x  = y / exmass
  sumup = 0.0D0
  sumdown = 0.0D0
  DO j = 1, iter
     sumup = sumup + x * DBLE(j**2) * &
          EXP(-x * DBLE(j**2))
     sumdown = sumdown + EXP(-x*DBLE(j**2))
  END DO
  delE = kB * Tr * (sumup / sumdown) ! in meV

END SUBROUTINE delta
