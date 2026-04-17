!>
!! \file gcd.f90
!! \brief return the greatest common divisor
!<

!>
!! \ingroup tbdftsnt
!! \brief Returns the greatest common divisor of two integer numbers
!!
!! Subroutine gcd returns the greatest common divisor (gcd) 
!! of two integers by using Euclidean algorithm.
!!
!! \param[in] n an integer
!! \param[in] m an integer
!! \param[out] d greatest common divisor of n and m
!!
!! \section References
!! -# R.Saito et al., Physical Properties of Carbon Nanotubes
!<
subroutine gcd( n, m, d )

  implicit none

  integer,intent(in) :: n, m
  integer,intent(out) :: d
  integer :: i, j, k, l

  i = iabs(n)
  j = iabs(m)
  if (j.gt.i) then
     k = j
     j = i
     i = k
  endif
  if (j.eq.0) then
     l = i
  else
     k = j
     do while (k.ne.0)
        k = mod(i,j)
        if (k.eq.0) then
           l = j
        else
           i = j
           j = k
        endif
     enddo
  endif

  d = l

end subroutine gcd
