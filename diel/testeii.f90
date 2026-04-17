program testeii

  implicit none

  integer :: nn, mm, ii, nkappa
  real(8) :: kappa, lk, dt, pdtlk
  
  real(8) :: Eii(40,0:30,4,80,2)
  real(8) :: dttheta(40,0:30,2)

  open(50,file="eii.dat",status="old",form="unformatted")
  read(50) Eii, dttheta
  close(50)

  nn = 6
  mm = 6
  ii = 1
  kappa = 2.2
 
  nkappa = int(10*(kappa-1)+1)
  WRITE(*,*) nkappa
  lk = Eii(nn,mm,ii,nkappa,2)
  write(*,*) lk
  dt = dttheta(nn,mm,1)
  write(*,*) dt
  pdtlk = ii**0.8 * (1/dt)**1.6 * (1/lk)**0.4
  print"(a,f8.5)", "pdtlk = ", pdtlk

end program testeii
