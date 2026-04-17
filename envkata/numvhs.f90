
integer function numvhs(n,m)

  implicit none

  integer,intent(in) :: n, m

  integer :: mu, line
  integer,parameter :: max = 1000

  line = 0
  do mu = 0, max
     select case(n-m)
     case(0)
        if( (n+m<4*mu).and.(4*mu<3*(n+m)) ) then
           line = line + 1
        end if
     case default
        if( (n<2*mu).and.(2*mu<2*n+m) ) then
           line = line + 1
        end if
     end select
  end do

  numvhs = line

end function numvhs
