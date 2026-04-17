program findoptkapp
! find an optimized kappa from experimental Eii value

  implicit none

  ! general variables
  integer :: j, k, dummy, opt
  integer, allocatable :: famNum(:) ! family from input to output
  integer, allocatable :: p(:) ! cutting line distance ratio
  integer :: nlinei ! number of line in input
  integer :: del ! for removing large differences
  ! number of line in output should be nlinei-del

  ! variables for reading input of experiments data
  ! and later write it to the output file
  integer, allocatable :: n(:), m(:), i(:)
  real(8), allocatable :: Ei(:)
  character(30) :: inputfile

  ! variables for writing output file
  real(8), allocatable :: dt(:), E(:), lk(:), optkapp(:) 
  character(30) :: outfile

  ! variables for reading sx-eXX.dat:
  integer :: nliner
  integer, allocatable :: nr(:), mr(:)
  real(8), allocatable :: dr(:), thetar(:), Er(:), lkr(:), kappr(:)
  real(8), allocatable :: diff(:)
  character(25) :: filetoread
  character(15), parameter :: root = "genkap/"
  
  print"(a)", "#Input file: "
  read*, inputfile
  ! inputfile format: n, m, i, E(exp)
  open(1,file=trim(inputfile),status="old")
  outfile = trim("opt-"//inputfile)
  open(3,file=outfile)
  nlinei = 0
  do
     read(1,*,end=10) dummy
     nlinei = nlinei+1
  end do
10 rewind(1)
  allocate(n(nlinei))
  allocate(m(nlinei))
  allocate(i(nlinei))
  allocate(Ei(nlinei))
  allocate(famNum(nlinei))
  allocate(dt(nlinei))
  allocate(E(nlinei))
  allocate(lk(nlinei))
  allocate(optkapp(nlinei))
  allocate(p(nlinei))

  del = 0 ! initializing number of data which have large energy difference
  do j = 1, nlinei
     read(1,*), n(j), m(j), i(j), Ei(j)
     famNum(j) = 2*n(j) + m(j)
     print*, mod(famNum(j),3)
     if (mod(famNum(j),3) == 0) then
        filetoread = trim(root)//"s0-e"//char(48+i(j))//char(48+i(j))//".dat"
        if (i(j) == 2 .or. i(j) == 3) then
           p(j) = 3
        else if (i(j) == 4 .or. i(j) == 5) then
           p(j) = 6
        else
           print*, "undefined p!"
        end if
     else if (mod(famNum(j),3) == 1) then
        filetoread = trim(root)//"s1-e"//char(48+i(j))//char(48+i(j))//".dat"
        if (i(j) == 1) then
           p(j) = 1
        else if (i(j) == 2) then
           p(j) = 2
        else if (i(j) == 3) then
           p(j) = 4
        else if (i(j) == 4) then
           p(j) = 5
        else
           print*, "undefined p!"
        end if
     else if (mod(famNum(j),3) == 2) then
        filetoread = trim(root)//"s2-e"//char(48+i(j))//char(48+i(j))//".dat"
        if (i(j) == 1) then
           p(j) = 1
        else if (i(j) == 2) then
           p(j) = 2
        else if (i(j) == 3) then
           p(j) = 4
        else if (i(j) == 4) then
           p(j) = 5
        else
           print*, "undefined p!"
        end if
     else
        print*, "incorrect family number!"
     end if
     
     ! file to read: n, m, dt, theta, E(cal), lk, kappa
     open(2,file=trim(filetoread),status="old")
     nliner = 0
     do
        read(2,*,end=20), dummy
        nliner = nliner+1
     end do
     print*, nliner
20   rewind(2)
     allocate(nr(nliner))
     allocate(mr(nliner))
     allocate(dr(nliner))
     allocate(thetar(nliner))
     allocate(Er(nliner))
     allocate(lkr(nliner))
     allocate(kappr(nliner))
     allocate(diff(nliner))
     do k = 1, nliner
        read(2,*), nr(k), mr(k), dr(k), thetar(k), Er(k), &
             & lkr(k), kappr(k)
        if (n(j) == nr(k) .and. m(j) == mr(k)) then
           diff(k) = abs(Ei(j)-Er(k))
        else
           diff(k) = 1.d0 ! good enough to be not considered
        end if
     end do
     opt = minloc(diff,dim=1)
     n(j) = nr(opt); m(j) = mr(opt); dt(j) = dr(opt)
     diff(j) = 1000*diff(opt)
     E(j) = Er(opt); lk(j) = lkr(opt)
     optkapp(j) = kappr(opt)
     if (abs(diff(j)<4.d0)) then
        write(3,"(2i6,3f10.3,f7.1,i6)"), n(j), m(j), dt(j), &
             & Ei(j), lk(j), optkapp(j), p(j)
     else
        del = del + 1
     end if
     deallocate(nr,mr,dr,thetar,Er,lkr,kappr,diff)
     close(2)
  end do
  
  close(1)
  close(3)
  
  deallocate(n,m,i,Ei,famNum,p)
  deallocate(dt,E,lk,optkapp)
  write(*,"(a,i3)"), "Number of lines in inputfile : ", nlinei  
  write(*,"(a,i3)"), "Number of lines in outputfile: ", nlinei-del
end program
