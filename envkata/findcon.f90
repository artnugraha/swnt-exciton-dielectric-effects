subroutine findcon(dataname,elkei,resultfunc)

implicit none

integer :: i, numline, locelkei
character :: dummy
character(100) :: dataname
real(8) :: elkei
real(8), allocatable :: sumbux(:), sumbuy(:), difference(:)
real(8), intent(out) :: resultfunc

open(11,file=dataname,status='old')
numline = 0
do
   read(11,*,end=21) dummy
   numline = numline+1
end do
21 rewind(11)

allocate(sumbux(numline))
allocate(sumbuy(numline))
allocate(difference(numline))

do i = 1,numline
   read(11,*) sumbux(i), sumbuy(i)
   difference(i) = abs(sumbux(i) - elkei)
end do
locelkei = minloc(difference,dim=1)
resultfunc = sumbuy(locelkei)
close(11)

deallocate(sumbux)
deallocate(sumbuy)
deallocate(difference)

end subroutine findcon
