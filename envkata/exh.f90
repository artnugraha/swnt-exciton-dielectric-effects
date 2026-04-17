!>
!! \file exh.f90
!! \brief Construct exciton hamiltonian and Calculate eigenvalue and eigenvectors
!<

!>
!! \brief Construct exciton hamiltonian and Calculate eigenvalue and eigenvectors
!!
!! \param[out] exval nth exciton energy (exval(1) correspods to 0th exciton energy)
!! \param[out] exvec eigen vector
!!
!! \note
!! -# LAPACK zhegv is used.
!! -# eigen vector is always calculated. See evflag parameter in this source code.
!!
!! \author Jie Jiang (jiang@chips.ncsu.edu)
!! \author Kentaro Sato (kentaro@flex.phys.tohoku.ac.jp)
!<
subroutine exh(otype,stype,n,nmkr,ec,ev,wc,wx,exval,exvec)

  use common,only: bound, nn, mm, nmrel

  implicit none

  !---------------------------------------------------------------------
  ! INPUT
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  character(len=2),intent(in) :: otype
  ! orbital type

  character(len=1),intent(in) :: stype
  ! spin type

  integer,intent(in) :: n
  ! dimension of Hamiltonian

  type(nmrel) :: nmkr

  complex(8),intent(in) :: ec(2,bound)
  ! quasi-particle energy (conduction band)

  complex(8),intent(in) :: ev(2,bound)
  ! quasi-particle energy (valence band)

  complex(8),intent(in) :: wc(2,2,bound,bound)
  ! direct Coulomb interaction

  complex(8),intent(in) :: wx(2,2,bound,bound)
  ! exchange Coulomb interaction

  !---------------------------------------------------------------------
  ! OUTPUT
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(8),intent(out) :: exval(bound)
  ! exciton energy

  complex(8),intent(out) :: exvec(bound,bound)
  ! exciton wave function coefficient

  !---------------------------------------------------------------------
  ! LOCAL
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer :: nk1, nk2
  complex(8) :: ecn, evn
  complex(8) :: wcn, wvn, wcn1, wcn2, wxn, wxn1, wxn2
  complex(8) :: h(bound,bound) ! Hamiltonian
  integer :: i1, j1, i2, j2, i, j
  complex(8),parameter :: c0 = ( 0.0d0, 0.0d0 ) ,c1 = ( 1.0d0, 0.0d0 )

  !---------------------------------------------------------------------
  ! variables for LAPACK zhegv
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer,parameter :: evflag = 1
  character(len=1) :: jobz
  character(len=1),parameter :: ceval ='N', cevec ='V', uplo = 'U'
  integer,parameter ::  itype = 1
  integer :: lda, ldb, lwork, info
  complex(8) :: a(n,n), b(n,n), work(2*n-1)
  real(8) :: w(n), rwork(3*n-2)

  lda = n
  ldb = n
  lwork = 2*n - 1

  nk1 = 0
  do i1 = 1, nmkr%nline
     do j1 = 1, nmkr%nk(i1)
        nk1 = nk1 + 1
        if( nk1.gt.n ) then
           write (*,100)
           stop
        endif

        nk2 = 0
        do i2 = 1, nmkr%nline
           do j2 = 1, nmkr%nk(i2)
              nk2 = nk2 + 1
              if( nk1.gt.n ) then
                 write (*,100)
                 stop
              endif

              ecn = 0.0d0
              evn = 0.0d0

              if( nk1.eq.nk2 ) then
                 select case(otype)
                 case('ge')
                    ! A1 exciton
                    !ecn = ec(1,nk1)
                    !evn = ev(1,nk1)

                    ecn = 0.5d0*(ec(1,nk1)+ec(2,nk1))
                    evn = 0.5d0*(ev(1,nk1)+ev(2,nk1))
                 case('go')
                    ! A2 exciton
                    !ecn = ec(1,nk1)
                    !evn = ev(1,nk1)

                    ecn = 0.5d0*(ec(1,nk1)+ec(2,nk1))
                    evn = 0.5d0*(ev(1,nk1)+ev(2,nk1))
                 case('k1')
                     ! E exciton
                    ecn = ec(1,nk1)
                    evn = ev(1,nk1)
                 case('k2')
                     ! E exciton
                    ecn = ec(1,nk1)
                    evn = ev(1,nk1)
                 end select
              endif

              select case(otype)
              case('ge')
                 ! A1 exciton
                 !wcn = wc(1,1,nk1,nk2)

                 wcn1 = 0.5d0*(wc(1,1,nk1,nk2)+wc(2,2,nk1,nk2))
                 wcn2 = 0.5d0*(wc(1,2,nk1,nk2)+wc(2,1,nk1,nk2))
                 wcn = wcn1 - wcn2
              case('go')
                 ! A2 exciton
                 !wcn = wc(1,1,nk1,nk2)

                 wcn1 = 0.5d0*(wc(1,1,nk1,nk2)+wc(2,2,nk1,nk2))
                 wcn2 = 0.5d0*(wc(1,2,nk1,nk2)+wc(2,1,nk1,nk2))
                 wcn = wcn1 + wcn2
              case('k1')
                 ! E exciton
                 wcn = wc(1,1,nk1,nk2)
              case('k2')
                 ! E exciton
                 wcn = wc(1,1,nk1,nk2)
              end select

              if( stype.eq.'s' ) then
                 select case(otype)
                 case ('ge')
                    ! A1 exciton
                    wxn = 0.0d0
                    ! From eq.(23), J. Jiang et al., PRB 75, 035407 (2007), 
                    ! singlet A1 state is always degenerated with triplet bonding state. 
                    ! Energy difference between A2 and A1 excitons
                    ! is determined by exchange term Kx from Coulomb interaction.
                    ! See also T. Ando, JPSJ 75, 024707 (2006).

                    !wxn1 = 0.5d0*(wx(1,1,nk1,nk2)+wx(2,2,nk1,nk2))
                    !wxn2 = 0.5d0*(wx(1,2,nk1,nk2)+wx(2,1,nk1,nk2))
                    !wxn = wxn1 - wxn2
                 case ('go')
                    ! A2 exciton
                    !wxn = 0.0d0

                    wxn1 = 0.5d0*(wx(1,1,nk1,nk2)+wx(2,2,nk1,nk2))
                    wxn2 = 0.5d0*(wx(1,2,nk1,nk2)+wx(2,1,nk1,nk2))
                    wxn = wxn1 + wxn2
                 case('k1')
                    ! E exciton
                    wxn = wx(1,1,nk1,nk2)
                 case('k2')
                    ! E exciton
                    wxn = wx(1,1,nk1,nk2)
                 end select
              elseif ( stype.eq.'t' ) then
                 wxn = 0.0d0
              endif

              h(nk1,nk2) = (ecn-evn) + 2.0d0*wxn - wcn
              ! left side of BS equation

           enddo
        enddo

     enddo
  enddo

  !=begin=lapack=diagonalization=!
  if(evflag.eq.0) then
     jobz = ceval
  else
     jobz = cevec
  endif

  do i = 1, n
     do j = 1, n
        a(j,i) = h(j,i)
     enddo
  enddo

  do i = 1, n
     do j = 1, i-1
        b(j,i) = c0
     enddo
     b(i,i) = c1
     do j = i+1, n
        b(j,i) = c0
     enddo
  enddo

  call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)

  if( info.ne.0 ) then
     ! diagonalization is failed
     write(*,200)
     stop
  endif

  ! matrix size of exval is bound by bound
  ! matrix size of w is n by n
  exval(1:n) = w(1:n)

  if( evflag.ne.0 ) then
     ! matrix size of exvec is bound by bound
     ! matrix size of a is n by n
     exvec(1:n,1:n) = a(1:n,1:n)
  endif
  !=end=lapack=diagonalization=!

100 format (" >> exh >> dimension problem")
200 format(" >> exh >> lapack >> zhegv >> info >> overflow >> ")

end subroutine exh
