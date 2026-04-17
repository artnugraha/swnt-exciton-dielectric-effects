
!*************************************************************************
!     construct exciton hamiltonian and calculate
!     eigenvalue and eigenvectors
!*************************************************************************
subroutine exh(otype,stype,n,nmkr,ec,ev,wc,wx,wc2,wx2,wc3,wx3,wcn2g,exval,exvec)

  use common,only:bound
  use common,only:nmrel
  use common,only:nn,mm
  use common,only : kexchange, kdirect, nc, nu

  use exkataif,only : num2str

  implicit none

  ! input
  integer,intent(in) :: n
  character(len=2),intent(in) :: otype
  character(len=1),intent(in) :: stype
  complex(8),intent(in) :: ec(2,bound), ev(2,bound), wc(2,2,bound,bound), wx(2,2,bound,bound)
  ! for short-range Coulomb parameter w1 and w2
  complex(8),intent(in) :: wc2(2,2,bound,bound), wx2(2,2,bound,bound), wc3(2,2,bound,bound), wx3(2,2,bound,bound)
  complex(8) :: wcn2g

  !output
  real(8) exval(bound)
  complex(8) exvec(bound,bound)
  !local
  integer nk1,nk2
  complex(8) ecn,evn
  complex(8) wcn,wvn,wcn1,wcn2,wxn,wxn1,wxn2
  complex(8) h(bound,bound)

  type(nmrel)::nmkr
  integer i1,j1,i2,j2,i,j

  integer*4 evflag
  parameter (evflag=1)
  complex(8) c0,c1
  parameter (c0=(0.0d0,0.0d0),c1=(1.0d0,0.0d0))
  !LAPACK variables for lapack zhegv start here
  character*1 jobz,ceval,cevec,uplo
  parameter (ceval='N',cevec='V',uplo='U')
  integer itype,lda,ldb,lwork,info
  parameter (itype=1)
  !lda,ldb
  complex(8) a(n,n)
  complex(8) b(n,n)
  real(8) w(n)
  !lwork
  complex(8) work(2*n-1)
  real(8) rwork(3*n-2)
  !LAPACK variables for lapack zhegv end here

  character(len=4) :: senbu
  character(len=100) :: fn500

  wcn2g = (0.0D0,0.0D0)

  kdirect(:) = (0.0D0,0.0D0)
  kexchange(:) = (0.0D0,0.0D0)

  senbu =  num2str(nn,mm)
  fn500 = './srcp/srcp-swnt'//senbu//'.dat'
  open(500,FILE=fn500)

  lda=n
  ldb=n
  lwork=2*n-1

  nk1=0
  do i1=1,nmkr%nline
     do j1=1,nmkr%nk(i1)
        nk1=nk1+1
        if (nk1.gt.n) then
           write (*,100)
           stop
        endif

        nk2=0
        do i2=1,nmkr%nline
           do j2=1,nmkr%nk(i2)
              nk2=nk2+1
              if (nk1.gt.n) then
                 write (*,100)
                 stop
              endif

              ecn=0.0d0
              evn=0.0d0

              if( nk1.eq.nk2 ) then
                 select case(otype)
                 case('ge')
                    ecn=ec(1,nk1)
                    evn=ev(1,nk1)

                    !ecn=0.5d0*(ec(1,nk1)+ec(2,nk1))
                    !evn=0.5d0*(ev(1,nk1)+ev(2,nk1))
                 case('go')
                    ecn=ec(1,nk1)
                    evn=ev(1,nk1)

                    !ecn=0.5d0*(ec(1,nk1)+ec(2,nk1))
                    !evn=0.5d0*(ev(1,nk1)+ev(2,nk1))
                 case('k1')
                    ecn=ec(1,nk1)
                    evn=ev(1,nk1)
                 case('k2')
                    ecn=ec(1,nk1)
                    evn=ev(1,nk1)
                 end select
              endif

              select case(otype)
              case('ge') ! A1 exciton
                 !wcn=wc(1,1,nk1,nk2)

                 wcn1=0.5d0*(wc(1,1,nk1,nk2)+wc(2,2,nk1,nk2))
                 wcn2=0.5d0*(wc(1,2,nk1,nk2)+wc(2,1,nk1,nk2))
                 wcn=wcn1-wcn2
              case('go') ! A2 exciton
                 !wcn=wc(1,1,nk1,nk2)

                 wcn1 = 0.5d0*(wc(1,1,nk1,nk2)+wc(2,2,nk1,nk2))
                 wcn2 = 0.5d0*(wc(1,2,nk1,nk2)+wc(2,1,nk1,nk2))
                 wcn = wcn1 + wcn2
              case('k1')
                 wcn=wc(1,1,nk1,nk2)
              case('k2')
                 wcn=wc(1,1,nk1,nk2)
              end select

              if (stype.eq.'s')then
                 select case(otype)
                 case ('ge') ! A1 exciton
                    wxn = 0.0d0 ! Eq.(23), J. Jiang et al., PRB 75, 035407 (2007).

                    !wxn1 = 0.5d0*(wx(1,1,nk1,nk2)+wx(2,2,nk1,nk2))
                    !wxn2 = 0.5d0*(wx(1,2,nk1,nk2)+wx(2,1,nk1,nk2))
                    !wxn = wxn1 - wxn2
                 case ('go') ! A2 exciton
                    !wxn=0.0d0

                    wxn1 = 0.5d0*(wx(1,1,nk1,nk2)+wx(2,2,nk1,nk2))
                    wxn2 = 0.5d0*(wx(1,2,nk1,nk2)+wx(2,1,nk1,nk2))
                    wxn = wxn1 + wxn2
                 case('k1')
                    wxn=wx(1,1,nk1,nk2)
                 case('k2')
                    wxn=wx(1,1,nk1,nk2)
                 end select
              elseif (stype.eq.'t') then
                 wxn=0.0d0
              endif

              h(nk1,nk2)=(ecn-evn)+2.0d0*wxn-wcn

              kexchange(1) = kexchange(1) + wx(1,1,nk1,nk2)
              kexchange(2) = kexchange(2) + wx(2,1,nk1,nk2)
              kexchange(3) = kexchange(3) + wx(1,2,nk1,nk2)
              kexchange(4) = kexchange(4) + wx(2,2,nk1,nk2)

              kdirect(1) = kdirect(1) + wc(1,1,nk1,nk2)
              kdirect(2) = kdirect(2) + wc(2,1,nk1,nk2)
              kdirect(3) = kdirect(3) + wc(1,2,nk1,nk2)
              kdirect(4) = kdirect(4) + wc(2,2,nk1,nk2)

!              write(500,FMT='(2I7,16F30.18)') nk1, nk2&
!                   &, wx(2,1,nk1,nk2), wx2(2,1,nk1,nk2), abs(wx(2,1,nk1,nk2)/wx2(2,1,nk1,nk2))&
!                   &, wx3(2,1,nk1,nk2), abs(wx3(2,1,nk1,nk2))&
!                   &, wx(1,2,nk1,nk2), wx2(1,2,nk1,nk2), abs(wx(1,2,nk1,nk2)/wx2(1,2,nk1,nk2))&
!                   &, wx3(1,2,nk1,nk2), abs(wx3(1,2,nk1,nk2))

!              write(500,FMT='(2I7,I15,16F25.10)') nk1, nk2, nc*nu&
!                   &, wx(1,1,nk1,nk2), wx(2,1,nk1,nk2), wx(1,2,nk1,nk2), wx(2,2,nk1,nk2)&
!                   &, wc(1,1,nk1,nk2), wc(2,1,nk1,nk2), wc(1,2,nk1,nk2), wc(2,2,nk1,nk2)

              wcn2g = wcn2g + wcn2
              write(500,FMT='(3F25.10,2I4)') wcn2, abs(wcn2), nn, mm

           enddo
        enddo

     enddo
  enddo

  close(500)

  !=begin=lapack=diagonalization=!
  if (evflag.eq.0) then
     jobz=ceval
  else
     jobz=cevec
  endif
  do i=1,n
     do j=1,n
        a(j,i)=h(j,i)
     enddo
  enddo
  do i=1,n
     do j=1,i-1
        b(j,i)=c0
     enddo
     b(i,i)=c1
     do j=i+1,n
        b(j,i)=c0
     enddo
  enddo
  !write (*,300)
  call zhegv (itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
  if (info.ne.0) then
     write(*,200)
     stop
  endif
  do i=1,n
     exval(i)=w(i)
  enddo
  if (evflag.ne.0) then
     do i=1,n
        do j=1,n
           exvec(j,i)=a(j,i)
        enddo
     enddo
  endif
  !=end=lapack=diagonalization=!

  return

100 format (" >> exh >> dimension problem")
200 format(" >> exh >> lapack >> zhegv >> info >> overflow >> ")
300 format("begin to diagonalize")

end subroutine exh
