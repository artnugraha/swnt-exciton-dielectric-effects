subroutine dgaussls(fname,fwhm)

! dgaussls - Gaussian Least Squares Fitting (double precision)
! last modified: 2009.05.13, by A.R.T Nugraha
! nugraha@flex.phys.tohoku.ac.jp
!
! ver.0.5, transform into a subroutine, 2009.11.19
! ver.0.4, fix the initial guess problem, 2009.05.13
! ver.0.3, add initial guess improvement, 2009.05.11
! ver.0.2, modify gaussls into dgaussls, 2009.05.06
! ver.0.1, original code, 2009.05.04
! 
! Purpose:
!     - find parameter \beta_1, \beta_2, \beta_3
!       !>
!          gaussian function of the form:
!          f(x) = A exp(-(x-x0)^2/(2\sigma^2))
!          \beta_1 --> A
!          \beta_2 --> x0
!          \beta_3 --> \sigma
!       <!
!     - calculate FWHM = 2*\sigma*sqrt(2*ln(2))
!
! Method:
!     - LeastSquares with GaussNewton algorithm
!       !>
!          data points:
!          (x1,y1), (x2,y2), ... (x_nlines, y_nlines)
!       !>
!          define residuals for i = 1,2, ..., nlines:
!          res_i = y_i - f(x_i;\beta), 
!          \sum {res_i} should be minimized
!       !>
!          linearized estimation:
!          starting with initial guess \beta^{0}, then
!          \beta^{j+1} = \beta^{j} + \delta{\beta}
!            (for j iterations)
!          and \delta{\beta} satisfies matrix equation:
!          (jacobian^T * jacobian) \delta{\beta} = &
!                & -(jacobian^T)*res
!       <!
!       if \delta{\beta} is small enough --> convergent

implicit none

! main variables
real(8), parameter :: eps0 = 1.0d-5
real(8), allocatable :: res(:), jacobi(:,:)
real(8), allocatable :: x(:), y(:)
real(8) :: jac(3,3), v(3), jacr(3), eps
real(8) :: beta1, beta2, beta3
real(8), intent(out) :: fwhm

! variables for beta3 guess:
real(8) :: xhalf, yhalf, xhalf1, xhalf2, yhalf1, yhalf2
integer :: lochalf1, lochalf2

! file operation, lapack, and iteration purpose
integer :: i, j, k, m, nlines, nrhs, info
integer :: ipiv(3)
character(50) :: fname

! read wavefunction data

open(1,file=fname,status='old')
nlines = 0
do
   read(1,*,end=10) eps ! temporary, use again later
   nlines = nlines+1
end do
10 rewind(1)
allocate(x(nlines))
allocate(y(nlines))
read(1,*) (x(i), y(i), i = 1,nlines)
close(1)

allocate(res(nlines))
allocate(jacobi(nlines,3))

! initialization (guess)
beta1 = maxval(y)-0.02; beta2 = x(maxloc(y,dim=1));

! guess for beta3 (i.e. sigma or standard deviation)
yhalf = 0.5 * maxval(y)
lochalf1 = maxloc(y,dim=1,mask=((y.lt.yhalf).and.(x.lt.beta2)))
lochalf2 = lochalf1 + 1

xhalf1 = x(lochalf1)
yhalf1 = y(lochalf1)
xhalf2 = x(lochalf2)
yhalf2 = y(lochalf2)
xhalf = xhalf1+(xhalf2-xhalf1)*((yhalf-yhalf1)/(yhalf2-yhalf1))

! beta3 = 2 * abs(beta2-xhalf)
beta3 = 0.68 * (2 * abs(beta2-xhalf))

! lapack parameters
m = 3 ! (m,n) matrix, m = n, lda = ldb = m
nrhs = 1

! matrix initialization
res = 0.; jacobi = 0.

! main algorithm
v = (/ beta1, beta2, beta3 /)

j = 0 ! initialize iteration
do    
   do k = 1, nlines
      res(k) = y(k) - fmain(beta1,beta2,beta3,x(k))
      jacobi(k,1) = -df1(beta2,beta3,x(k))
      jacobi(k,2) = -df2(beta1,beta2,beta3,x(k))
      jacobi(k,3) = -df3(beta1,beta2,beta3,x(k))
   end do
   jac = matmul(transpose(jacobi),jacobi)
   jacr = -matmul(transpose(jacobi),res)
   
   ! solving matrix equation
   ! (jacobian^T * jacobian) \delta{\beta} = -(jacobian^T)*res
   ! analogous to AX = B, use lapack subroutine dgetrf and dgetrs
   ! (or sgetrf/sgetrs for single precision)
   
   call dgetrf(m,m,jac,m,ipiv,info)
   call dgetrs('N',m,nrhs,jac,m,ipiv,jacr,m,info)
   ! jacr now becomes \delta{\beta}
   ! check convergence for loop exit
   eps = sqrt(jacr(1)**2 + jacr(2)**2 + jacr(3)**2)
   !write(*,'(1X,I4,4F9.5)') j, v(1), v(2), v(3), eps   
   if (eps < eps0) exit

   ! continue, \beta = \beta + \delta{\beta}
   v = v + jacr
   beta1 = v(1); beta2 = v(2); beta3 = v(3)
   j = j + 1 ! iteration

end do

fwhm = 2.0d0 * beta3 * sqrt(2.0d0*log(2.))

deallocate(x)
deallocate(y)

contains
! jacobian, partial derivatives of f(x) respect to \beta
  function fmain(b1,b2,b3,x)
    real(8) :: fmain
    real(8), intent(in) :: b1, b2, b3, x
    fmain = b1*exp(-(x-b2)**2/(2*b3**2))
  end function fmain

  function df1(b2,b3,x)
    real(8) :: df1
    real(8), intent(in) :: b2, b3, x
    df1 = exp(-(x-b2)**2/(2*b3**2))
  end function df1
 
  function df2(b1,b2,b3,x)
    real(8) ::  df2
    real(8), intent(in) :: b1, b2, b3, x
    df2 = (b1*(x-b2)/b3**2) * exp(-(x-b2)**2/(2*b3**2))
  end function df2

  function df3(b1,b2,b3,x)
    real(8) ::  df3
    real(8), intent(in) :: b1, b2, b3, x
    df3 = (b1*(x-b2)**2/b3**3) * exp(-(x-b2)**2/(2*b3**2))
  end function df3

end subroutine dgaussls
