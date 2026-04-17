subroutine wfcoecb(nmwf1,nmwf2,nmwf3,wfcoe3,wfcoe4)
use common,only:wfmesh
use common,only:wfqp
use common,only:nc,nt
use common,only:nab,rab,hab,oab
implicit none
!input
type(wfmesh)::nmwf1,nmwf2,nmwf3
!output
type(wfqp)::wfcoe3,wfcoe4
!local
!tbtube
integer*4 nb,nd,nu,ns
parameter (nb=4,nd=3,nu=64,ns=2)
integer*4 nm
parameter (nm=nb*ns)
integer*4 evflag
parameter (evflag=1)
real*8 ener(nm)
complex*16 ampl(nm,nm)
!
integer mubt,muup,kbt,kup
real*8 ev,ec
complex*16 cv(2),cc(2)
integer muk,k
! for allocate
integer sev,sec,scv,scc
logical cond
!

mubt=nmwf2%mubt
muup=nmwf2%muup
kbt=nmwf2%kbt
kup=nmwf2%kup

allocate (wfcoe3%ev(mubt:muup,kbt:kup),stat=sev)
allocate (wfcoe3%ec(mubt:muup,kbt:kup),stat=sec)
allocate (wfcoe3%cv(mubt:muup,kbt:kup,2),stat=scv)
allocate (wfcoe3%cc(mubt:muup,kbt:kup,2),stat=scc)

cond=sev.eq.0.and.sec.eq.0.and.scv.eq.0.and.scc.eq.0

if (.not.cond ) then
write (*,20)
endif


do muk=nmwf2%mubt,nmwf2%muup,1
do k=nmwf2%kbt,nmwf2%kup,1
call tbtube (evflag,muk,k,nc,nt,nab,rab,hab,oab,ener,ampl)
ev=ener(4)
ec=ener(5)
cv(1)=ampl(2,4)
cv(2)=ampl(6,4)
cc(1)=ampl(2,5)
cc(2)=ampl(6,5)
wfcoe3%ev(muk,k)=ev
wfcoe3%ec(muk,k)=ec
wfcoe3%cv(muk,k,:)=cv(:)
wfcoe3%cc(muk,k,:)=cc(:)
enddo
enddo

mubt=nmwf3%mubt
muup=nmwf3%muup
kbt=nmwf3%kbt
kup=nmwf3%kup

allocate (wfcoe4%ev(mubt:muup,kbt:kup),stat=sev)
allocate (wfcoe4%ec(mubt:muup,kbt:kup),stat=sec)
allocate (wfcoe4%cv(mubt:muup,kbt:kup,2),stat=scv)
allocate (wfcoe4%cc(mubt:muup,kbt:kup,2),stat=scc)

cond=sev.eq.0.and.sec.eq.0.and.scv.eq.0.and.scc.eq.0

if (.not.cond ) then
write (*,30)
endif


do muk=nmwf3%mubt,nmwf3%muup,1
do k=nmwf3%kbt,nmwf3%kup,1
call tbtube (evflag,muk,k,nc,nt,nab,rab,hab,oab,ener,ampl)
ev=ener(4)
ec=ener(5)
cv(1)=ampl(2,4)
cv(2)=ampl(6,4)
cc(1)=ampl(2,5)
cc(2)=ampl(6,5)
wfcoe4%ev(muk,k)=ev
wfcoe4%ec(muk,k)=ec
wfcoe4%cv(muk,k,:)=cv(:)
wfcoe4%cc(muk,k,:)=cc(:)
enddo
enddo


20 format ("allocate problem for wfcoe3")
30 format ("allocate problem for wfcoe4")


end subroutine wfcoecb
