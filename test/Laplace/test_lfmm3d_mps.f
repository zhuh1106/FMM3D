       implicit none
       integer :: ns, nt, nmpole
       integer :: i,j,k,ntest,nd,idim,ier,ijk,ilen,isuccess,iper
       integer :: ifcharge,ifdipole,ifpgh,ifpghtarg
       integer :: ipass(18),len1,ntests,isum,lused,lw,n1,nlege
       integer :: npts, ns1, ntarg, ntm, ntot
       integer, allocatable :: mterms(:), impole(:)
       integer :: interms,lca
       integer, allocatable :: tmp_vec(:)
        
       double precision :: eps, err, hkrand, dnorm, done, h, pi
       double precision :: rscale, sc, shift, thresh
       double precision, allocatable :: source(:,:), targ(:,:)
       double precision, allocatable :: cmpole(:,:), dc(:,:)
       double precision, allocatable :: wlege(:), rmpole(:)
        
       double precision, allocatable :: charge(:,:)
       double precision, allocatable :: dipvec(:,:,:)
       double precision, allocatable :: pot(:,:),pot2(:,:),pottarg(:,:)
       double precision, allocatable :: grad(:,:,:),gradtarg(:,:,:)
       double precision, allocatable :: hess(:,:,:),hesstarg(:,:,:)
       double complex, allocatable :: mpole(:), local(:)
       double complex, allocatable :: ilocal(:,:)


cccccc        
cccccc        
cccccc        
cccccc        
cccccc temporary variable cccccccccccccccccccccccccccccccccccccccccccccccc 
cccccc   
c
cc     tree variables
c
       integer idivflag, ndiv, nboxes, nlevels
       integer nlmax, nlmin, ifunif
       integer *8 ipointer(8),ltree
       integer ipointer4(8),ltree4
       integer, allocatable :: itree(:)
       integer, allocatable :: isrcse(:,:),isrc(:)
       integer itargtree, ntargtree
       double precision :: targtree(3)
       double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc     temporary sorted arrays
c
       integer :: lmpole, mt
       integer, allocatable :: mtermssort(:), impolesort(:)
       double precision, allocatable :: cmpolesort(:,:)
       double precision, allocatable :: rmpolesort(:)
       double complex, allocatable :: mpolesort(:)
       double complex, allocatable :: localsort(:)

c
cc     temporary fmm arrays
c
       integer, allocatable :: nterms(:)
       integer *8, allocatable :: iaddr(:,:)
       double precision, allocatable :: scales(:)
       double precision, allocatable :: rmlexp(:)

       integer lmptemp,nmax
       integer *8 lmptot
       double precision, allocatable :: mptemp(:),mptemp2(:)


c
cc     other temporary variables
c
       integer iert,ifprint,ilev,l
       double precision time1,time2,omp_get_wtime
       integer, allocatable :: laddr(:,:)

c      temporary variables for fmps main subroutine
c      temp variables       
       integer, allocatable :: nlist1(:),list1(:,:)
       integer, allocatable :: nlist2(:),list2(:,:)
       integer, allocatable :: nlist3(:),list3(:,:)
       integer, allocatable :: nlist4(:),list4(:,:)

       double precision timeinfo(6)
       integer isep

       integer nuall,ndall,nnall,nsall,neall,nwall
       integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
       integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
       integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

       integer, allocatable :: uall(:,:),dall(:,:),nall(:,:)
       integer, allocatable :: sall(:,:),eall(:,:),wall(:,:)
       integer, allocatable :: u1234(:,:),d5678(:,:)
       integer, allocatable :: n1256(:,:),s3478(:,:)
       integer, allocatable :: e1357(:,:),w2468(:,:)
       integer, allocatable :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
       integer, allocatable :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
       integer, allocatable :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
       integer, allocatable :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

c      temp variables
       integer ii,jj,kk,ll,m,iloc
       integer ibox,jbox,npts0
       integer nchild

       integer istart,iend
       integer jstart,jend

       double precision d

c      PW variables
       integer nlams, nthmax, nphmax, nmaxt
       double precision, allocatable :: carray(:,:)
       double precision, allocatable :: rdplus(:,:,:)
       double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
       double precision, allocatable :: rdmsq3(:,:,:)
       double precision, allocatable :: rlams(:),whts(:)

       double precision, allocatable :: rlsc(:,:,:)
       integer, allocatable :: nfourier(:), nphysical(:)
       integer nexptot, nexptotp
       double complex, allocatable :: xshift(:,:)
       double complex, allocatable :: yshift(:,:)
       double precision, allocatable :: zshift(:,:)

       double complex, allocatable :: fexpe(:),fexpo(:),fexpback(:)
       double complex, allocatable :: mexp(:,:,:,:)
       double complex, allocatable :: mexpf1(:,:,:),mexpf2(:,:,:)
       double complex, allocatable :: mexpp1(:,:,:),mexpp2(:,:,:)
       double complex, allocatable :: mexppall(:,:,:,:)

       double complex, allocatable :: tmp(:,:,:,:)
       double precision, allocatable :: mptmp(:,:)

       double precision rtmp

       integer lw7, lused7

       integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
       integer nn
       double precision, allocatable :: rscpow(:)

c      list 3 variables
       double complex, allocatable :: iboxlexp(:,:,:)
       double precision, allocatable :: iboxsubcenters(:,:,:)
       integer, allocatable :: iboxsrcind(:,:)
       integer, allocatable :: iboxisort(:,:)
       integer, allocatable :: iboxisort_tmp(:,:)
       integer, allocatable :: iboxfl(:,:,:)
c      end of list 3 variables

c      list 4 variables
       integer cntlist4
       integer, allocatable :: list4ct(:),ilist4(:)
       double complex, allocatable :: gboxmexp(:,:,:)
       double complex, allocatable :: gboxwexp(:,:,:,:,:)
       double complex, allocatable :: pgboxwexp(:,:,:,:)
       double precision, allocatable :: gboxsubcenters(:,:,:)
       integer, allocatable :: gboxind(:,:)
       integer, allocatable :: gboxisort(:,:)
       integer, allocatable :: gboxisort_tmp(:,:)
       integer, allocatable :: gboxfl(:,:,:)
c      end of list 4 variables


       integer *8 bigint

       integer nthd,ithd
       integer omp_get_max_threads,omp_get_thread_num   

cccccc        
cccccc        
cccccc        
cccccc       
cccccc set up cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cccccc         
       lca = 4*50
       done = 1
       pi = 4*atan(done)

c      initialize printing routine
       call prini(6,13)
      
       nd = 1
      
      
       n1 = 4
       ns = n1**3
       nmpole = ns
      
       call prinf('ns = *', ns, 1)
       call prinf('nc = *', nmpole, 1)
        
       nt = 19
      
       allocate(source(3,ns),targ(3,nt), cmpole(3,nmpole))
       allocate(charge(nd,ns),dipvec(nd,3,ns))
       allocate(pot(nd,ns), pot2(nd,ns))
       allocate(grad(nd,3,ns))
       allocate(hess(nd,6,ns))

       allocate(pottarg(nd,nt))
       allocate(gradtarg(nd,3,nt))
       allocate(hesstarg(nd,6,nt))
       eps = 0.5d-9
      
       write(*,*) "=========================================="
       write(*,*) "Testing suite for lfmm3d_mps"
       write(*,'(a,e12.5)') "Requested precision = ",eps
      
c      generate sources uniformly in the unit cube 
       h = 1.0d0/(n1+1)
       ijk = 0
       do i = 1,n1
         do j = 1,n1
           do k = 1,n1
             ijk = ijk + 1
             source(1,ijk) = h*i
             source(2,ijk) = h*j
             source(3,ijk) = h*k
           enddo
         enddo
       enddo
        
        
       dnorm = 0
       do i=1,ns
         do idim=1,nd
           charge(idim,i) = hkrand(0)
           dnorm = dnorm + abs(charge(idim,i))**2

           dipvec(idim,1,i) = hkrand(0)
           dipvec(idim,2,i) = hkrand(0)
           dipvec(idim,3,i) = hkrand(0)
      
           pot(idim,i) = 0
           grad(idim,1,i) = 0
           grad(idim,2,i) = 0
           grad(idim,3,i) = 0
         enddo
       enddo
      
       dnorm = sqrt(dnorm)
       do i=1,ns
         do idim = 1,nd
           charge(idim,i) = charge(idim,i)/dnorm
         end do
       end do
        
      
       shift = h/1000
       call prin2('shift = *', shift, 1)
       do i = 1,ns
         cmpole(1,i) = source(1,i) + shift
         cmpole(2,i) = source(2,i)
         cmpole(3,i) = source(3,i)
       end do
      
c      now form a multipole expansion at each center
       allocate(mterms(nmpole), impole(nmpole))
      
       ntm = 7
       ntot = 0
       do i = 1,nmpole
         mterms(i) = ntm
         ntot = ntot + (mterms(i)+1)*(2*mterms(i)+1)
       end do
      
       allocate(mpole(nd*ntot))
      
       impole(1) = 1
       do i = 1,nmpole-1
         ilen = (mterms(i)+1)*(2*mterms(i)+1)
         impole(i+1) = impole(i) + nd*ilen
       end do
      
       
       nlege = 300
       lw = 5*(nlege+1)**2
       allocate( wlege(lw) )
      
       call prinf('before ylgndrfwini, lw = *', lw, 1)
       call ylgndrfwini(nlege, wlege, lw, lused)
       call prinf('after ylgndrfwini, lused = *', lused, 1)
      
       call zinitialize(nd*ntot*2, mpole)
       
       ns1 = 1
       rscale = 1
       sc = shift
       if (sc .lt. 1) rscale = sc
      
       call prin2('rscale = *', rscale, 1)
       
       allocate( rmpole(nmpole) )

       do i = 1,nmpole
         rmpole(i) = rscale
         call l3dformmpc(nd, rscale, source(1,i), charge(1,i),
     1         ns1, cmpole(1,i), mterms(i), mpole(impole(i)),
     2         wlege, nlege)
       enddo

c      do the direct calculation
       thresh = 1.0d-15
       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ntarg = 0
       ifpghtarg = 0
       ier = 0
       call lfmm3d(nd, eps, ns, source, ifcharge,
     1      charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess,
     2      ntarg, targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier)
       
       call prin2('via fmm, potential = *', pot, 10)
      
       
       allocate(local(nd*ntot))
       local = 0

cccccc        
cccccc        
cccccc        
cccccc 
cccccc mps call cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
cccccc      
c      now test source to source, charge,with potentials
       print *
       print *
       print *
       write(6,*) 'testing multipoles to locals'
       write(6,*) 'input: multipole expansions'
       write(6,*) 'output: local expansions'
       write(6,*) 
       write(6,*) 
    !    call lfmm3dmps0(nd, eps,
    !  1      nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier)
c
c
c      ifprint is an internal information printing flag. 
c      Suppressed if ifprint=0.
c      Prints timing breakdown and other things if ifprint=1.
c      
call cpu_time(time1)
C$     time1=omp_get_wtime()      
       ifprint=1

       print *, 'ndiv still needs to be optimized'
       ndiv = 1

       if(ifprint.ge.1) print *, "ndiv = ", ndiv
c
cc     figure out tree structure
c
c
cc     set tree flags
c 
       idivflag = 0
       nlmax = 51
       nlevels = 0
       nboxes = 0
       ltree = 0
       nlmin = 0
       iper = 0
       ifunif = 0
c
cc
c
       iert = 0
       ntargtree = 0
       targtree(1) = 0
       targtree(2) = 0
       targtree(3) = 0

c
cc     memory management code for contructing level restricted tree
       call ptstreemem(cmpole,nmpole,targtree,ntargtree,idivflag,ndiv,
     1      nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree4)
c      to avoid integer 8 mex crash...     
       ltree = int(ltree4, 8)
      
       allocate(itree(ltree))
       allocate(boxsize(0:nlevels))
       allocate(treecenters(3,nboxes))

c       Call tree code
c      input order of ifunif&iper is not consistent with pts_tree_build 
c      subroutine
       call ptstreebuild(cmpole,nmpole,targtree,ntargtree,idivflag,
     1      ndiv,nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree4,itree,
     2      ipointer4,treecenters,boxsize)
       ipointer = int(ipointer4, 8)
      

       allocate(isrcse(2,nboxes))
       allocate(isrc(nmpole))

       call ptstreesort(nmpole,cmpole,itree,ltree4,nboxes,nlevels,
     1      ipointer4,treecenters,isrc,isrcse)
      
c
cc     end of tree build
c
       if(ifprint.ge.1) print *, "ltree/1e9", ltree/1.0d9
       if(ifprint.ge.1) print *, "nlevels=", nlevels


       allocate(nterms(0:nlevels))
c      Compute length of expansions at each level
       nmax = 0
       do i=0,nlevels
         call l3dterms(eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
       enddo
c       
c      Multipole and local expansions will be held in workspace
c      in locations pointed to by array iaddr(2,nboxes).
c
c      iaddr is pointer to iaddr array, itself contained in workspace.
c      imptemp is pointer for single expansion (dimensioned by nmax)
c
c      ... allocate iaddr and temporary arrays
c

       allocate(iaddr(2,nboxes))
       lmptemp = (nmax+1)*(2*nmax+1)*2*nd
       allocate(mptemp(lmptemp),mptemp2(lmptemp))

c
cc     reorder pole centers
c
       allocate(cmpolesort(3, nmpole))
       allocate(rmpolesort(nmpole))
       allocate(impolesort(nmpole))
       allocate(mtermssort(nmpole))

       lmpole = 0
       do i = 1,nmpole
         lmpole = lmpole + (mterms(i)+1)*(2*mterms(i)+1)
       enddo
       lmpole = nd * lmpole

       allocate(mpolesort(lmpole))

       call dreorderf(3, nmpole, cmpole, cmpolesort, isrc)
       call dreorderf(1, nmpole, rmpole, rmpolesort, isrc)
       call ireorderf(1, nmpole, mterms, mtermssort, isrc)

       impolesort(1) = 1
       do i=1,nmpole
         mt = mtermssort(i)
         ilen = (mt+1)*(2*mt+1)

         ijk = 1
         do j = 1, ilen
           do l = 1,nd
             mpolesort(impolesort(i)+ijk-1)=mpole(impole(isrc(i))+ijk-1)
             ijk = ijk + 1
           enddo
         enddo

         if(i.lt.nmpole) impolesort(i+1) = impolesort(i) + nd*ilen
       enddo

c
c      allocate memory need by multipole, local expansions at all
c      levels
c      irmlexp is pointer for workspace need by various fmm routines,
c
       call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
       if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9
       allocate(rmlexp(lmptot),stat=iert)
       if(iert.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         ier = 4
         return
       endif

       allocate(localsort(lmpole))
       localsort = 0

c      Memory allocation is complete.
c      scaling factor for multipole and local expansions at all levels
c
       allocate(scales(0:nlevels))
       do ilev = 0,nlevels
         scales(ilev) = boxsize(ilev)
       enddo

       call cpu_time(time2)
C$     time2=omp_get_wtime()      

       if(ifprint.ge.1) 
     1   call prin2('time before fmm main=*',time2-time1,1)
c     Call main fmm routine

       call cpu_time(time1)
       ier = 0
       allocate(laddr(2,0:nlevels))
       laddr = reshape(itree(ipointer(1):2*(nlevels+1)),(/2,nlevels+1/))
      
c     Call main fmm routine

       call cpu_time(time1)
       ier = 0
C$     time1=omp_get_wtime()
       call lfmm3dmain_mps0(nd,eps,
     $            nmpole,cmpolesort,rmpolesort,mtermssort,mpolesort,
     $            impolesort,localsort,
     $            iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $            itree,ltree,ipointer,ndiv,nlevels,
     $            nboxes,iper,boxsize,treecenters,isrcse,
     $            scales,laddr,nterms,ier)
       if(ier.ne.0) return
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(6) = time2-time1
       if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
       d = 0
       do i = 1,6
         d = d + timeinfo(i)
       enddo
       if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
       if(ier.ne.0) return
       

       do i = 1,nmpole
         mt = mtermssort(i)
         ilen = (mt+1)*(2*mt+1)

         ijk = 1
         do j = 1,ilen
           do l = 1,nd
             local(impole(isrc(i))+ijk-1)=localsort(impolesort(i)+ijk-1)
             ijk = ijk+1
           enddo
         enddo
       enddo       

cccccc        
cccccc        
cccccc        
cccccc 
cccccc post process cccccccccccccccccccccccccccccccccccccccccccccccccccccc   
cccccc 
       interms = (mterms(1)+1)*(2*mterms(1)+1)
       allocate(ilocal(mterms(1)+1,-mterms(1):mterms(1)))
       ilocal = 0
       allocate(tmp_vec(interms))
       call zinitialize(nd*nmpole, pot2)
       npts = 1
       do i = 1,nmpole
         call l3dtaevalp(nd, rmpole(i),
     1        cmpole(1,i), local(impole(i)),
     2        mterms(i), source(1,i), npts, pot2(1,i),
     3        wlege, nlege)
       enddo
       print *, size(local)
       
       call prin2('from lfmm3d_mps, potential = *', pot2, 10)
      
       err = 0
       dnorm = 0
       do j = 1,nmpole
         do i = 1,nd
           err = err + abs(pot(i,j)-pot2(i,j))**2
           dnorm = dnorm + abs(pot(i,j))**2
           pot2(i,j) = pot2(i,j) - pot(i,j)
         end do
       end do
       
       err = sqrt(err/dnorm)
       call prin2('l2 rel err=*',err,1)
      
       open(unit=33,file='print_testres.txt',access='append')
       isuccess = 0
       ntest = 1
       if(err.lt.eps) isuccess = 1
      
       write(33,'(a,i1,a,i1,a)') 'Successfully completed ',
     1    isuccess,' out of ',ntest,' tests in helm3d_mps testing suite'
       close(33)
      
       stop
       end
      
       subroutine zinitialize(len, zs)
         implicit double precision (a-h,o-z)
         double precision :: zs(len)
      
         do i = 1,len
           zs(i) = 0
         end do
         return
       end subroutine zinitialize
