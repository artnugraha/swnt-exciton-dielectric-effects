!>
!! \file cutoff.f90
!! \brief cutoff function for self-energy
!<

!>
!! \brief cutoff function for self-energy
!!
!! \param[in] x is energy (conduction or valence band)
!! \param[out] cutoff
!!
!! \section References
!! -# Eq.(2.12) H. Sakai et al., J.Phys.Soc.Jpn. 72, 1698 (2003).
!!  - http://dx.doi.org/10.1143/JPSJ.72.1698 
!! -# Eq.(3.8), T.Ando, J.Phys.Soc.Jpn. 75, 024707 (2006).
!!  - http://dx.doi.org/10.1143/JPSJ.75.024707
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
function cutoff(x)

  use common,only : dt
  use exkataparameter,only : epc, alp

  implicit none

  real(8),intent(in) :: x
  real(8) :: cutoff
  real(8) :: ep

  ep = epc * sqrt(3.0d0) * 3.0d0 / dt
  cutoff = ep**alp / ( abs(x)**alp + ep**alp )

end function cutoff
