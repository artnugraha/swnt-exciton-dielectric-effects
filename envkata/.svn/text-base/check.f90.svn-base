
subroutine checkmu(type,mucon,muval,index,nmkr)

  use common,only:nmrel

  implicit none

  ! input
  character(len=2),intent(in) :: type
  integer,intent(in) :: mucon(2), muval(2), index
  type(nmrel),intent(in) :: nmkr

  ! local
  integer muconcom(2),muvalcom(2)
  logical cond
  integer i

  i=index
  muconcom(:)=(/nmkr%muc(i),-nmkr%muv(i)/)
  muvalcom(:)=(/nmkr%muv(i),-nmkr%muc(i)/)
  if (type.eq.'ga') then
     cond=mucon(1).ne.muconcom(1).or.mucon(2).ne.muconcom(2).or.&
          muval(1).ne.muvalcom(1).or.muval(2).ne.muvalcom(2)
  else
     cond=mucon(1).ne.muconcom(1).or.&
          muval(1).ne.muvalcom(1)
  endif
  if (cond) then
     write (*,50)
     stop
  endif

50 format (">>ini>>con or val wave vector (mu) is not correct")

end subroutine checkmu



subroutine checkk(type,kcon,kval,index1,index2,nmkr)

  use common,only:nmrel

  implicit none

  ! input
  character(len=2),intent(in) :: type
  integer,intent(in) :: kcon(2), kval(2), index1, index2
  type(nmrel),intent(in) :: nmkr

  !local
  integer kconcom(2),kvalcom(2)
  logical cond
  integer i,j

  i=index1
  j=index2

  kconcom(:)=(/nmkr%kc(i,j),-nmkr%kv(i,j)/)
  kvalcom(:)=(/nmkr%kv(i,j),-nmkr%kc(i,j)/)

  if (type.eq.'ga') then
     cond=kcon(1).ne.kconcom(1).or.kcon(2).ne.kconcom(2).or.&
          kval(1).ne.kvalcom(1).or.kval(2).ne.kvalcom(2)
  else
     cond=kcon(1).ne.kconcom(1).or.&
          kval(1).ne.kvalcom(1)
  endif
  if (cond) then
     write (*,50)
     stop
  endif

50 format (">>ini>>con or val wave vector (k) is not correct")

end subroutine checkk



subroutine checkrange(kc,kcbt,kcup,kt,ktbt,ktup,kvapex,sub)

  use common,only : nn, mm

  ! input
  integer,intent(in) :: kc, kcbt, kcup, kt, ktbt, ktup
  character(len=2),intent(in) :: kvapex
  character(len=3),intent(in) :: sub

  ! local
  logical :: cond

  cond = (kc.ge.kcbt).and.(kc.le.kcup).and.(kt.ge.ktbt).and.(kt.le.ktup)
  if(.not.cond) then
     write(*,50) kvapex, sub
     write(*,*) (kc.ge.kcbt), (kc.le.kcup), (kt.ge.ktbt), (kt.le.ktup)
     write(*,*) kc, kcbt, kcup, kt, ktbt, ktup
     write(*,*) nn, mm
     stop
  end if

50 format (A2,A3,"wavevector range is wrong")

end subroutine checkrange
