
module exkataif
  != interfaces for exkatampi.f90

  implicit none

  interface

     subroutine cutline(center,nmke)
       use common,only : nmmesh
       ! input
       integer,intent(in) :: center
       ! output
       type(nmmesh),intent(out) :: nmke
     end subroutine cutline
     

     subroutine qline(nmk1,nmq)
       use common,only : nmmesh, qmesh
       ! input
       type(nmmesh),intent(in) :: nmk1
       ! output
       type(qmesh),intent(out) :: nmq
     end subroutine qline
  

     subroutine ftdia(muc,kc,nmkr,nmk1,nmk2,nmq,&
          nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,con,con11,&
          con12,con21,vq,vqc11,vqc12,vqc21)
       use common,only : nmrel, nmmesh,qmesh, wfqp, wfmesh
       ! input
       integer,intent(in) :: muc, kc
       type(nmrel),intent(in) :: nmkr
       type(nmmesh),intent(in) :: nmk1, nmk2
       type(qmesh),intent(in) :: nmq
       ! output
       type(wfmesh),intent(out) :: nmwf2, nmwf3
       type(wfqp),intent(out) :: wfcoe1, wfcoe2, wfcoe3
       real(8),pointer :: con(:,:), con11(:,:), con12(:,:), con21(:,:)
       complex(8),pointer :: vq(:,:,:,:), vqc11(:,:,:,:), vqc12(:,:,:,:), vqc21(:,:,:,:)
     end subroutine ftdia
 

     subroutine wflinecbi(nmk1,nmk2,nmwf1,nmwf2)
       use common,only : nmmesh, wfmesh
       ! input
       type(nmmesh),intent(in) :: nmk1, nmk2
       ! output
       type(wfmesh),intent(out) :: nmwf1, nmwf2
     end subroutine wflinecbi


     subroutine wfcoecbi(nmwf1,nmwf2,wfcoe1,wfcoe2)
       use common,only : wfmesh, wfqp
       ! input
       type(wfmesh),intent(in) :: nmwf1, nmwf2
       ! output
       type(wfqp),intent(out) :: wfcoe1, wfcoe2
     end subroutine wfcoecbi
     

     subroutine mukrg(nline1,mu1,bt1,up1,nline2,mu2,bt2,up2,nline,mu,bt,up)
       use common,only : mu2d, mumx
       ! input
       integer,intent(in) :: nline1, mu1(mu2d), bt1(mu2d), up1(mu2d)
       integer,intent(in) :: nline2, mu2(mu2d), bt2(mu2d), up2(mu2d)
       ! output
       integer,intent(out) :: nline, mu(mumx), bt(mumx), up(mumx)
     end subroutine mukrg


     subroutine wflinecb(nmk1,nmk2,nmwf1,nmwf2,nmwf3)
       use common,only : nmmesh, wfmesh
       ! input
       type(nmmesh),intent(in) :: nmk1, nmk2
       ! output
       type(wfmesh),intent(out) :: nmwf1, nmwf2, nmwf3
     end subroutine wflinecb
  

     subroutine wfcoecb(nmwf1,nmwf2,nmwf3,wfcoe3,wfcoe4)
       use common,only : wfmesh, wfqp
       ! input
       type(wfmesh) :: nmwf1, nmwf2, nmwf3
       ! output
       type(wfqp) :: wfcoe3, wfcoe4
     end subroutine wfcoecb


     subroutine diacb(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,wfcoe4,vq2,con)
       use common,only : nmmesh, wfqp, wfmesh
       ! input
       integer,intent(in) :: muq, q
       type(nmmesh),intent(in) :: nmk1, nmk2
       type(wfmesh),intent(in) :: nmwf2, nmwf3
       type(wfqp),intent(in) :: wfcoe1, wfcoe2, wfcoe3, wfcoe4
       complex(8),intent(in) :: vq2(2,2)
       ! output
       real(8),intent(out) :: con
     end subroutine diacb
     

     subroutine selfcb(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,wfcoe4,selfenergy)
       use common,only : nmmesh, wfmesh, wfqp
       ! input
       integer,intent(in) :: muq, q
       type(nmmesh),intent(in) :: nmk1, nmk2
       type(wfmesh),intent(in) :: nmwf2, nmwf3
       type(wfqp),intent(in) :: wfcoe1, wfcoe2, wfcoe3, wfcoe4
       ! output
       real(8),intent(out) :: selfenergy
     end subroutine selfcb
     

     subroutine wfline(nmk1,nmk2,nmwf1,nmwf2,nmwf3)
       use common,only : nmmesh, wfmesh
       ! input
       type(nmmesh),intent(in) :: nmk1, nmk2
       type(wfmesh),intent(in) :: nmwf1
       ! nmwf1 is not used in subroutine wfline. This variable must be intent(out).
       ! output
       type(wfmesh),intent(out) :: nmwf2, nmwf3
     end subroutine wfline

     
     subroutine wfcoe(nmwf1,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3)
       use common,only : wfmesh, wfqp
       ! input
       type(wfmesh),intent(in) :: nmwf1, nmwf2, nmwf3
       ! output
       type(wfqp),intent(out) :: wfcoe1, wfcoe2, wfcoe3
     end subroutine wfcoe


     subroutine dia(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,vq2,con)
       use common,only : nmmesh, wfqp, wfmesh
       ! input
       integer,intent(in) :: muq, q
       complex(8),intent(in) :: vq2(2,2)
       type(nmmesh),intent(in) :: nmk1, nmk2
       type(wfmesh),intent(in) :: nmwf2, nmwf3
       type(wfqp),intent(in) :: wfcoe1, wfcoe2, wfcoe3
       ! output
       real(8),intent(out) :: con
     end subroutine dia
     

     function cutoff(x)
       ! input
       real(8) :: x
       ! output
       real(8) :: cutoff
     end function cutoff
     
 
     subroutine self(muq,q,nmk1,nmk2,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,selfenergy)
       use common,only : nmmesh, wfmesh, wfqp
       ! input
       integer,intent(in) :: muq, q
       type(nmmesh),intent(in) :: nmk1, nmk2
       type(wfmesh),intent(in) :: nmwf2, nmwf3
       type(wfqp),intent(in) :: wfcoe1, wfcoe2, wfcoe3
       ! output
       real(8),intent(out) :: selfenergy
     end subroutine self
     

     subroutine iniga(muc,kc,nmkr,nmq,con,con11,con12,con21&
          &,vq,vqc11,vqc12,vqc21,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3&
          &,ec,ev,ecs,evs,wc,wx,kcv,fcv)
       use common,only : bound, nmrel,qmesh, wfqp, wfmesh
       ! input
       integer,intent(in) :: muc, kc
       type(nmrel),intent(in) :: nmkr
       type(qmesh),intent(in) :: nmq
       real(8),pointer :: con(:,:), con11(:,:), con12(:,:), con21(:,:)
       complex(8),pointer :: vq(:,:,:,:), vqc11(:,:,:,:), vqc12(:,:,:,:), vqc21(:,:,:,:)
       type(wfmesh),intent(in) :: nmwf2, nmwf3
       type(wfqp),intent(in) :: wfcoe1, wfcoe2, wfcoe3
       ! output
       complex(8),intent(out) :: ec(2,bound), ev(2,bound), wc(2,2,bound,bound), wx(2,2,bound,bound)
       real(8),intent(out) :: ecs(2,bound), evs(2,bound)
       integer,intent(out) :: kcv(bound,2,2,2)
       complex(8),intent(out) :: fcv(bound,2,2,2)
     end subroutine iniga


     subroutine ftrans(mu,q,type,vq)
       use common,only : ns
       ! input
       integer(4),intent(in) :: mu, q
       integer(4),intent(in) :: type
       ! output
       complex(8),intent(out) :: vq(ns,ns)
     end subroutine ftrans


     subroutine checkmu(type,mucon,muval,index,nmkr)
       use common,only:nmrel
       ! input
       character(len=2),intent(in) :: type
       integer,intent(in) :: mucon(2), muval(2), index
       type(nmrel),intent(in) :: nmkr
     end subroutine checkmu


     subroutine checkk(type,kcon,kval,index1,index2,nmkr)
       use common,only : nmrel
       ! input
       character(len=2),intent(in) :: type
       integer,intent(in) :: kcon(2), kval(2), index1, index2
       type(nmrel),intent(in) :: nmkr
     end subroutine checkk
     

     subroutine checkrange(kc,kcbt,kcup,kt,ktbt,ktup,kvapex,sub)
       ! input
       integer,intent(in) :: kc, kcbt, kcup, kt, ktbt, ktup
       character(len=2),intent(in) :: kvapex
       character(len=3),intent(in) :: sub
     end subroutine checkrange

     
     subroutine qpencal(mukc,kc,mukv,kv,cc1,cv1,vq1,con,nmq,nmwf2,nmwf3,wfcoe1,wfcoe2,wfcoe3,qev,qec)
       use common,only : wfqp, wfmesh, qmesh
       ! input
       integer,intent(in) :: mukc, kc, mukv, kv
       complex(8),intent(in) :: cc1(2), cv1(2)
       complex(8),pointer :: vq1(:,:,:,:)
       real(8),pointer :: con(:,:)
       type(qmesh),intent(in) :: nmq
       type(wfmesh),intent(in) :: nmwf2, nmwf3
       type(wfqp),intent(in) :: wfcoe1, wfcoe2, wfcoe3
       ! output
       complex(8),intent(out) :: qev(2,2), qec(2,2)
     end subroutine qpencal


     subroutine exh(otype,stype,n,nmkr,ec,ev,wc,wx,exval,exvec)
       use common,only : bound, nmrel
       ! input
       character(len=2),intent(in) :: otype
       character(len=1),intent(in) :: stype
       integer,intent(in) :: n
       type(nmrel),intent(in) :: nmkr
       complex(8),intent(in) :: ec(2,bound), ev(2,bound), wc(2,2,bound,bound), wx(2,2,bound,bound)
       ! output
       real(8) :: exval(bound)
       complex(8) :: exvec(bound,bound)
     end subroutine exh


     function eg3(x,y)
       use exkataparameter,only : EGresult
       real(8),intent(in) :: x, y
       type(EGresult) :: eg3
     end function eg3


     subroutine findkc(type,nmk1,nmk2,ao,cline,vline,nmkc)
       use common,only : nmmesh, nmcen
       ! input
       character(len=2),intent(in) :: type
       character(len=3),intent(in) :: ao
       integer,intent(in) :: cline, vline
       type(nmmesh),intent(in) :: nmk1, nmk2
       ! output
       type(nmcen),intent(out) :: nmkc
     end subroutine findkc
     
     
     subroutine findkr(type,nmk1,nmk2,muc,kc,ao,cline,vline,nmkr)
       use common,only : nmmesh, nmrel
       ! input
       integer,intent(in) :: muc, kc, cline, vline
       character(len=2),intent(in) :: type
       character(len=3),intent(in) :: ao
       type(nmmesh),intent(in) :: nmk1, nmk2
       ! output
       type(nmrel),intent(out) :: nmkr
     end subroutine findkr
     

     character(4) function num2str(nn,mm)
       integer,intent(in) :: nn, mm
     end function num2str

     integer function numvhs(n,m)
       integer,intent(in) :: n, m
     end function numvhs


     subroutine inik(muc,kc,nmkr,ec,ev,wc,wx)
       use common,only : bound, nmrel
       integer :: muc, kc
       type(nmrel),intent(in) :: nmkr
       character(len=2) :: type
       complex(8) :: ec(2,bound), ev(2,bound), wc(2,2,bound,bound), wx(2,2,bound,bound)
     end subroutine inik


     subroutine SelectKappa(tubetype,cline,dt,kappa)
       character(len=2),intent(in) :: tubetype
       integer,intent(in) :: cline
       real(8),intent(in) :: dt
       real(8),intent(out) :: kappa
     end subroutine SelectKappa
     
     subroutine gcd( n, m, d )
       integer,intent(in) :: n, m
       integer,intent(out) :: d
     end subroutine gcd
  end interface

end module exkataif
