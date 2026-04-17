!>
!! \file common.f90
!! \brief global variables for exkatampi system
!<

!>
!! \brief global variables for exkatampi system
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
module common

  !> max array size for muarry1 etc. in cutline.90
  integer,parameter :: mu2d = 50

  integer,parameter :: mumx=102, murmx=50, krmx=1000, bound=2000

  integer :: nmukr
  real(8) v0

  !> the energy cost to place two electrons on a single site
  !! See INPUT PARAMETER section in exkatampi.f90
  real(8) :: u0

  ! a static dielectric constant describing the effects of electrons in core states
  ! See INPUT PARAMETER section in exkatampi.f90
  real(8) :: kapp


  !> cutting line and energy band information around K or K' point
  type :: nmmesh
     !> number of array (or number of cutting line) in the region.
     !> See cutline.f90.
     integer :: nline

     !> cutting line number (mu) in the region. See cutline.f90.
     integer :: mu(mu2d)

     !> k_vector at starting point on the cutting line mu
     integer :: bt(mu2d)

     !> k_vector at ending point on the cutting line mu
     integer :: up(mu2d)

     !> k_vector of vHSs point on the cutting line mu
     integer :: vhsk(mu2d)

     !> minimum energy on the cuttingline mu
     real(8) :: vhse(mu2d)

     !> minimum conduction band energy
     real(8) :: vhsec(mu2d)

     !> minimum valence band energy
     real(8) :: vhsev(mu2d)
  end type nmmesh


  !> qmesh
  type :: qmesh
     integer :: nline
     integer :: mu(mumx)
     integer :: bt, up, point
  end type qmesh


  !> bounds of wave vector around K or K' point for calculation. See also wflinecbi.f90.
  type :: wfmesh
     !> number of array (or number of cutting line)
     integer :: nline

     !> minimum number of nmmesh\%mu(:)
     integer :: mubt

     !> maximum number of nmmesh\%mu(:)
     integer :: muup

     !> minimum number of nmmesh\%bt(:)
     integer :: kbt

     !> minimum number of nmmesh\%up(:)
     integer :: kup

     !> number of points between kbt and kup
     integer :: point
  end type wfmesh


  !> wave function coefficient and energy around K or K point
  type :: wfqp
     !> pi (valence) band energy
     real(8),pointer :: ev(:,:)

     !> anti-pi (conduction) band energy
     real(8),pointer :: ec(:,:)

     !> wave function coefficient of pi (valence) band
     complex(8),pointer :: cv(:,:,:)

     !> wave function coefficient of anti-pi (conduction) band
     complex(8),pointer :: cc(:,:,:)
  end type wfqp


  !> nmcen
  type :: nmcen
     integer :: nline
     integer :: mu(mumx), bt(mumx), up(mumx)
  end type nmcen


  !> nmrel
  type :: nmrel
     integer :: nline
     integer,dimension(murmx) :: mu, muc, muv, nk
     integer,dimension(murmx,krmx) :: k, kc, kv
  end type nmrel


  !> number of type of neighbor site (A or B atom)
  integer,parameter :: ns = 2

  integer,parameter :: nu = 64000, nu1 = 1000
  integer :: nabf(ns,ns)
  real(8) :: xyzc(3,ns)
  real(8) :: xyzn(3,nu,ns,ns)
  real(8) :: rabc(2,ns)
  real(8) :: rabn(2,nu,ns,ns)


  ! t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nc,nt
  ! parameters for subroutine tubpar, tbini and tbtube.
  ! See also above source code.
  
  !> nanotube unit cell length [nm]
  real(8) :: tt
 
  !> nanotube diameter [nm]
  real(8) :: dt

  !> nanotube chiral angle [radian]
  real(8) :: theta
  
  !> translational distance between A & B atoms [t]
  real(8) :: tab

  !> circumferential distance between A & B atoms [radian]
  real(8) :: phiab

  !> dab  radial distance between A & B atoms [d]
  real(8) :: dab

  !> translational component of graphene unit vector a1 [t]
  real(8) :: a1t

  !> circumferential component of a1 vector [radian]
  real(8) :: a1phi

  !> translational component of graphene unit vector a2 [t]
  real(8) :: a2t

  !> circumferential component of a2 vector [radian]
  real(8) :: a2phi

  !> number of cutting lines in circumferential direction
  integer :: nc

  !> number of mesh points in translational direction
  integer :: nt
    
  !> number of atom pairs of a given type
  integer :: nab(ns,ns)

  !> interatomic distances for each atom pair
  real(8) :: rab(3,64,ns,ns)

  !> hamiltonian matrix blocks
  real(8) :: hab(4,4,64,ns,ns)

  !> overlap matrix blocks
  real(8) :: oab(4,4,64,ns,ns)

  integer :: unmax
  real(8) :: rabc1(2,ns)
  real(8) :: rabn1(2,-nu1:nu1,-nu1:nu1,ns,ns)

  !> apex for 2D BZ. ex. K1 or K2
  character(len=2) :: kcapex

  !> apex for 2D BZ. ex. K1 or K2
  character(len=2) :: kvapex

  !> chiral index m of (n,m)
  integer :: nn

  !> chiral index m of (n,m)
  integer :: mm

  !> reciplocal lattice vector K1 for (n,m) nanotube
  real(8) :: kbase1(2)

  !> reciplocal lattice vector K2 for (n,m) nanotube
  real(8) :: kbase2(2)

end module common
