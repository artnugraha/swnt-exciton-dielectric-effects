!>
!! \file exkataparameter.f90
!! \brief parameters for exkatampi system
!<

!>
!! \brief parameters for exkatampi system
!!
!! \section References
!! -# R.Saito et al., Physical Properties of Carbon Nanotubes.
!!
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
module exkataparameter

  implicit none

  !> PI
  real(8),parameter :: PI = 3.141592653589793D0

  !> sqrt(3)
  real(8),parameter :: SR3 = 1.732050807568877D0

  !> sqrt(2)
  real(8),parameter :: SR2 = 1.414213562373095D0

  !> imaginary unit
  complex(8),parameter :: I_U = ( 0.0D0, 1.0D0 )

  !> C-C bond lenght [angstrom]
  real(8),parameter :: ACC = 1.42D0

  !> latice constant [angstrom]
  real(8),parameter :: AUC = acc * sr3

  !> coupling parameter Hpp.
  !! See Table 2.1, R. Saito et al., Physical Properties of Carbon Nanotubes (MGM s617)
  real(8),parameter :: TPPP = -3.033D0

  !> coupling parameter Spp.
  !! See Table 2.1, R. Saito et al., Physical Properties of Carbon Nanotubes (MGM s617)
  real(8),parameter :: SPPP = 0.129D0

  !> Parameter for cutoff function.
  !! It is essentially independent of the parameter alp as long as 
  !! the cutoff function decays smoothly but rapidly enough.
  !! See also cutoff.f90 and references. 
  integer,parameter :: alp = 4

  !> Parameter for cutoff function. 
  !! epc is approximatly 3 gamma (gamma is transfer integral between nearest neighbor atomss).
  !! See also cutoff.f90 and references.
  real(8),parameter :: epc = 10.0d0

  !> Bohr radius [nm]
  real(8),parameter :: bohr = 0.5291772108d-01
  !> cutoff parameter for ETB. See tbcoord.f90.
  real(8) :: cutlength = 0.8d3*bohr

  !> K points in Brillouin zone of graphene
  real(8),parameter :: K_1(2) = (/  2.0D0*PI/SR3,  2.0D0*PI/3.0D0 /)
  real(8),parameter :: K_2(2) = (/         0.0D0,  4.0D0*PI/3.0D0 /)
  real(8),parameter :: K_3(2) = (/ -2.0D0*PI/SR3,  2.0D0*PI/3.0D0 /)
  real(8),parameter :: K_4(2) = (/ -2.0D0*PI/SR3, -2.0D0*PI/3.0D0 /)
  real(8),parameter :: K_5(2) = (/         0.0D0, -4.0D0*PI/3.0D0 /)
  real(8),parameter :: K_6(2) = (/  2.0D0*PI/SR3, -2.0D0*PI/3.0D0 /)

  !> M points in Brillouin zone of graphene
  real(8),parameter :: M_1(2) = (/  2.0D0*PI/SR3, 0.0D0 /)
  real(8),parameter :: M_2(2) = (/        PI/SR3,    PI /)
  real(8),parameter :: M_3(2) = (/       -PI/SR3,    PI /)
  real(8),parameter :: M_4(2) = (/ -2.0D0*PI/SR3, 0.0D0 /)
  real(8),parameter :: M_5(2) = (/       -PI/SR3,   -PI /)
  real(8),parameter :: M_6(2) = (/        PI/SR3,   -PI /)

  !> object composition for simple tight binding
  type :: EGresult
     !> energy. see eq.(2.27) in MGM s617
     !! w1:: valence band
     real(8) :: w1

     !> energy. see eq.(2.27) in MGM s617
     !! w2:: conduction band
     real(8) :: w2

     !> coefficients of wave function
     !! ca:: valence band
     complex(8) :: ca(2)

     !> coefficients of wave function
     !! cb:: conduction band
     complex(8) :: cb(2)
  end type EGresult

end module exkataparameter
