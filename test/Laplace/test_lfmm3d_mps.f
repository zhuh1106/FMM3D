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
       integer ipointer(8),ltree
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
       integer, allocatable :: iaddr(:,:)
       integer, allocatable :: iaddrtmp(:,:)
       double precision, allocatable :: scales(:)
       double precision, allocatable :: rmlexp(:)

       integer lmptemp,nmax
       integer lmptot
       integer lmptottmp
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
       double complex, allocatable :: mexpf12(:,:,:,:)
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


       integer bigint

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
      
       n1 = 7
      !  n1 = 9
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
           charge(idim,i) = real(i, 8)/ns
          !  charge(idim,i) = hkrand(0)
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
        
      
       shift = h/10
       call prin2('shift = *', shift, 1)
       do i = 1,ns
         cmpole(1,i) = source(1,i) + shift
         cmpole(2,i) = source(2,i) - shift
         cmpole(3,i) = source(3,i)
       end do
      
c      now form a multipole expansion at each center
       allocate(mterms(nmpole), impole(nmpole))
      
       ntm = 10
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

      !  print *, "wlege size : ", size(wlege)
      !  open(1, file = 'mps_data.dat')  
      !  do i=1,nd*ntot
      !   write(1,*) imag(mpole(i))
      !  enddo
      !  close(1)

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
      !  print *, "isrc : ", isrc
      !  print *, "isrcse : ", isrcse


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
      !  print *, "lmpole : ", lmpole
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
       allocate(iaddrtmp(2,nboxes))
       call mpalloc0(nd,itree(ipointer(1)),iaddrtmp,nlevels,lmptottmp,
     1      nterms,nboxes)
       print *, "lmptottmp = ", lmptottmp/1.0d9
       print *, "nterms = ", nterms
       lmptot = int(lmptottmp,8)
       iaddr = int(iaddrtmp,8)
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
       laddr = reshape(itree(ipointer(1):(ipointer(1)+2*(nlevels+1)-1)),
     1                 (/2,nlevels+1/))
      
c     Call main fmm_mps routine
      !  call cpu_time(time1)
       ier = 0
! C$     time1=omp_get_wtime()
!        call lfmm3dmain_mps0(nd,eps,
!      $            nmpole,cmpolesort,rmpolesort,mtermssort,mpolesort,
!      $            impolesort,localsort,
!      $            iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
!      $            itree,ltree,ipointer,ndiv,nlevels,
!      $            nboxes,iper,boxsize,treecenters,isrcse,
!      $            scales,laddr,nterms,ier)
!        if(ier.ne.0) return
!        call cpu_time(time2)
! C$     time2=omp_get_wtime()
       nthd = 1
C$     nthd=omp_get_max_threads()
       pi = 4.0d0*atan(1.0d0)
       thresh = 2.0d0**(-51)*boxsize(0)
       ifprint=1
c      initialize various tree lists
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 27
       isep = 1
       call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1      treecenters,itree(ipointer(3)),itree(ipointer(4)),
     2      itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2      itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
       print *, "mnlist4 = ", mnlist4
       allocate(list1(mnlist1,nboxes),nlist1(nboxes))
       allocate(list2(mnlist2,nboxes),nlist2(nboxes))
       allocate(list3(mnlist3,nboxes),nlist3(nboxes))
       allocate(list4(mnlist4,nboxes),nlist4(nboxes))

       print *, "list4: ", list4(1,1:4)

       call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1      treecenters,itree(ipointer(3)),itree(ipointer(4)),
     2      itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3      itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4      mnlist2,list2,nlist3,mnlist3,list3,
     4      nlist4,mnlist4,list4)
       
      open(1, file = 'mps_data8.dat')
      do i=1,577
        write(1,*) (list4(1,i))
      enddo
      close(1)  
      print *, "list4: ", list4(1,1:4)
       
c      Initialize routines for plane wave mp loc translation
       if(isep.eq.1) then
         if(eps.ge.0.5d-2) nlams = 12
         if(eps.lt.0.5d-2.and.eps.ge.0.5d-3) nlams = 12
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 20
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 29
         if(eps.lt.0.5d-9) nlams = 37
       endif
       if(isep.eq.2) then
         if(eps.ge.0.5d-3) nlams = 9
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 15
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 22
         if(eps.lt.0.5d-9) nlams = 29
       endif
       allocate(rlams(nlams),whts(nlams))
       allocate(nphysical(nlams),nfourier(nlams))
       nmax = 0
       do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
       enddo
       allocate(rscpow(0:nmax))
       allocate(carray(4*nmax+1,4*nmax+1))
       allocate(dc(0:4*nmax,0:4*nmax))
       allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
       allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
       allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
       allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
       allocate(rlsc(0:nmax,0:nmax,nlams))
c      generate rotation matrices and carray (rdminus was not initialized to 0?)
       call getpwrotmat(nmax,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)
c      generate rlams and weights (these are the nodes
c      and weights for the lambda integral)
       call vwts(rlams,whts,nlams)
c      generate the number of fourier modes required to represent the
c      moment function in fourier space
       call numthetahalf(nfourier,nlams)
c      generate the number of fourier modes in physical space
c      required for the exponential representation
       call numthetafour(nphysical,nlams)
c      Generate powers of lambda for the exponential basis
       call rlscini(rlsc,nlams,rlams,nmax)
c      Compute total number of plane waves
       nexptotp = 0
       nexptot = 0
       nthmax = 0
       nphmax = 0
       nn = 0
       do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
         nn = nn + nphysical(i)*nfourier(i)
       enddo
       print *, "nn size: ", nn
       allocate(fexpe(nn),fexpo(nn),fexpback(nn))
       allocate(tmp(nd,0:nmax,-nmax:nmax,nthd))
       allocate(mptmp(lmptemp,nthd))
       allocate(xshift(-5:5,nexptotp))
       allocate(yshift(-5:5,nexptotp))
       allocate(zshift(5,nexptotp))
       allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd),
     1          mexpp1(nd,nexptotp,nthd))
       allocate(mexpp2(nd,nexptotp,nthd),mexppall(nd,nexptotp,16,nthd))
       allocate(mexpf12(nd,nexptot,nthd,2))
cc     NOTE: there can be some memory savings here
       bigint = 0
       bigint = nboxes
       bigint = bigint*6
       bigint = bigint*nexptotp*nd
       if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9
       allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
       if(iert.ne.0) then
         print *, "Cannot allocate pw expansion workspace"
         print *, "bigint=", bigint
         ier = 8
         return
       endif
       allocate(list4ct(nboxes))
       allocate(ilist4(nboxes))
       do i=1,nboxes
         list4ct(i)=0
         ilist4(i)=0
       enddo
       cntlist4=0
c      Precompute table for shifting exponential coefficients in 
c      physical domain
       call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)
c      Precompute table of exponentials for mapping from
c      fourier to physical domain
       call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
       if(ifprint.ge.1) 
     1   call prin2('end of generating plane wave info*',i,0)
c      set timing to 0
       do i=1,6
         timeinfo(i)=0
       enddo
c      set all multipole and local expansions to zero
       do ilev = 0,nlevels
         do ibox=laddr(1,ilev),laddr(2,ilev)
           call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
           call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo       
       enddo
c      initialize legendre function evaluation routines
       nlege = 100
       lw7 = 5*(nlege+1)**2
       call ylgndrfwini(nlege,wlege,lw7,lused7)
cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 0 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
c      count number of boxes are in list4
       lca = 4*nmax
       if(ifprint.ge.1)
     $   call prinf('=== STEP 0 list4===*',i,0)
       call cpu_time(time1)
C$     time1=omp_get_wtime()
       do ilev=1,nlevels-1
         do ibox=laddr(1,ilev),laddr(2,ilev)
           if(nlist3(ibox).gt.0) then
             cntlist4=cntlist4+1
             list4ct(ibox)=cntlist4
             ilist4(cntlist4)=ibox
           endif
         enddo
       enddo
       if(ifprint.ge.1) print *,"nboxes:",nboxes,"cntlist4:",cntlist4
       allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
       allocate(gboxmexp(nd*(nterms(nlevels)+1)*
     1                   (2*nterms(nlevels)+1),8,cntlist4))
       allocate(gboxsubcenters(3,8,nthd))
       allocate(gboxfl(2,8,nthd))
       nmaxt = 0
       do ibox=1,nboxes
         if(list4ct(ibox).gt.0) then
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = iend-istart+1
           if(npts.gt.nmaxt) nmaxt = npts
         endif
       enddo
       print *, "nmaxt = ", nmaxt
       allocate(gboxind(nmaxt,nthd))
       allocate(gboxisort(nmaxt,nthd))
       allocate(gboxisort_tmp(nmaxt,nthd))
       allocate(gboxwexp(nd,nexptotp,6,8,nthd))
c      note gboxmexp is an array not scalar
       pgboxwexp=0d0
       gboxmexp=0d0
c      form mexp for all list4 type box at first ghost box center
       do ilev=1,nlevels-1
         rscpow(0) = 1.0d0/boxsize(ilev+1)
         rtmp = scales(ilev+1)/boxsize(ilev+1)
         do i=1,nterms(ilev+1)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc use n1 = 7 or 9 for list4 type box
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(list4ct(ibox).gt.0) then
              print *, "list4ct is processed"
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1
              if(npts.gt.0) then
                call subdividebox(cmpolesort(1,istart),npts,
     1               treecenters(1,ibox),boxsize(ilev+1),
     2               gboxind(1,ithd),gboxfl(1,1,ithd),
     3               gboxsubcenters(1,1,ithd))
                do i=istart,iend
                  gboxisort_tmp(i-istart+1,ithd) = i
                enddo
                call ireorderf(1,npts,gboxisort_tmp(1,ithd),
     1               gboxisort(1,ithd),gboxind(1,ithd))
                ! print *, "gboxisort_tmp = ", gboxisort_tmp
                ! print *, "gboxisort = ", gboxisort
                ! print *, "gboxind = ", gboxind
                ! print *, "gboxfl = ", gboxfl
                do i=1,8
                  if(gboxfl(1,i,ithd).gt.0) then
                    jstart=gboxfl(1,i,ithd)
                    jend=gboxfl(2,i,ithd)
                    npts0=jend-jstart+1
                    jbox=list4ct(ibox)
                    do j = jstart,jend
                      k = gboxisort(j,ithd)
                      call l3dmpmp(nd,rmpolesort(k),cmpolesort(1,k),
     1                     mpolesort(impolesort(k)),mtermssort(k),
     2                     scales(ilev+1),gboxsubcenters(1,i,ithd),
     3                     gboxmexp(1,i,jbox),nterms(ilev+1),dc,lca)
                    enddo
                    call l3dmpmp(nd,scales(ilev+1),
     1                   gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox),
     2                  nterms(ilev+1),scales(ilev),treecenters(1,ibox),
     3                   rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox),
     1                   rscpow,tmp(1,0,-nmax,ithd))
cc                process up down for current box
                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd),
     1                   nterms(ilev+1),nlams,
     2                   nfourier,nexptot,mexpf1(1,1,ithd),
     3                   mexpf2(1,1,ithd),rlsc)
                    call ftophys(nd,mexpf1(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,1,i,ithd),
     3                   fexpe,fexpo)
                    call ftophys(nd,mexpf2(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,2,i,ithd),
     3                   fexpe,fexpo)
                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd),
     1                   gboxwexp(1,1,2,i,ithd),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2),
     3                   xshift,yshift,zshift)
cc                process north-south for current box
                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),
     1                   mptmp(1,ithd),rdminus)
                    call mpoletoexp(nd,mptmp(1,ithd),
     1                   nterms(ilev+1),nlams,
     2                   nfourier,nexptot,mexpf1(1,1,ithd),
     3                   mexpf2(1,1,ithd),rlsc)
                    call ftophys(nd,mexpf1(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,3,i,ithd),
     3                   fexpe,fexpo)
                    call ftophys(nd,mexpf2(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,4,i,ithd),
     3                   fexpe,fexpo)
                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd),
     1                   gboxwexp(1,1,4,i,ithd),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4),
     3                   xshift,yshift,zshift)
cc                process east-west for current box
                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),
     1                   mptmp(1,ithd),rdplus)
                    call mpoletoexp(nd,mptmp(1,ithd),
     1                   nterms(ilev+1),nlams,
     2                   nfourier,nexptot,mexpf1(1,1,ithd),
     3                   mexpf2(1,1,ithd),rlsc)
                    call ftophys(nd,mexpf1(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,5,i,ithd),
     3                   fexpe,fexpo)
                    call ftophys(nd,mexpf2(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,6,i,ithd),
     3                   fexpe,fexpo)                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd),
     1                   gboxwexp(1,1,6,i,ithd),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6),
     3                   xshift,yshift,zshift)
                  endif
                enddo
              endif
            endif
         enddo
        !  if (ilev.eq.1) then
        !  open(1, file = 'mps_data12.dat')
        !   do i=1,6
        !     do j=1,1627
        !       write(1,*) real(pgboxwexp(1,j,1,i))
        !     enddo
        !   enddo
        ! endif
       enddo
      !  open(1, file = 'mps_data.dat')  
      !  do i=1,lmpole
      !   write(1,*) real(mpolesort(i))
      !  enddo
      !  close(1)
      !  open(1, file = 'mps_data2.dat')  
      !  do i=1,lmpole
      !   write(1,*) imag(mpolesort(i))
      !  enddo
      !  close(1)
      !  open(1, file = 'mps_data3.dat')  
      !  do i=1,nmpole
      !   write(1,*) rmpolesort(i)
      !  end do  
      !  close(1)
      !  open(1, file = 'mps_data4.dat')  
      !  do i=1,nmpole
      !   do j=1,3
      !     write(1,*) cmpolesort(j, i)
      !   enddo
      !  end do  
      !  close(1)
       open(1, file = 'mps_data12.dat')
       do i=1,6
         do j=1,1627
           write(1,*) real(pgboxwexp(1,j,1,i))
         enddo
       enddo

       
      !   open(1, file = 'mps_data.dat')  
      !  do i=1,1
      !   do j=-nmax,nmax
      !     do k=0,nmax
      !      write(1,*) imag(tmp(1,k,j,i))
      !     enddo
      !   enddo
      !  enddo
      !  close(1)

       deallocate(gboxfl,gboxsubcenters,gboxwexp)
       deallocate(gboxind)
       deallocate(gboxisort_tmp)
       deallocate(gboxisort)
       deallocate(gboxmexp)
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       if(ifprint.ge.1) print *,"mexp list4 time:",time2-time1
       timeinfo(3)=time2-time1
c      end of count number of boxes are in list4
cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 1 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
       if(ifprint .ge. 1)
     $   call prinf('=== STEP 1 (shift mp) ====*',i,0)
         call cpu_time(time1)
C$       time1=omp_get_wtime()
c       step 1, shift incoming multipole expansion to the center
c       of each leaf-node box
       print *, "size of rmlexp : ", size(rmlexp)
       do ilev=2,nlevels
         do ibox=laddr(1,ilev),laddr(2,ilev)
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = iend-istart+1
           nchild = itree(ipointer(4)+ibox-1)
           if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
             do i = istart,iend
               ! rmlexp is defined as real but is treated as complex inside l3dmpmp
               call l3dmpmp(nd,rmpolesort(i),cmpolesort(1,i),
     1              mpolesort(impolesort(i)),mtermssort(i),
     2              scales(ilev),treecenters(1,ibox),
     3              rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
             enddo
           endif
         enddo         
       enddo
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       
      !  open(1, file = 'mps_data.dat')  
      !  do i=1,size(rmlexp)
      !     write(1,*) rmlexp(i)
      !  end do  
      !  close(1)    
      
cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 2 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc       
       if(ifprint .ge. 1)
     $   call prinf('=== STEP 2 (merge mp) ====*',i,0)
       call cpu_time(time1)
C$     time1=omp_get_wtime()
       do ilev=nlevels-1,0,-1
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = isrcse(1,jbox)
                  iend = isrcse(2,jbox)
                  npts = iend-istart+1
                  if(npts.gt.0) then
                     call l3dmpmp(nd,scales(ilev+1),
     1               treecenters(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),scales(ilev),treecenters(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
                  endif
               endif
            enddo
         enddo    
       enddo
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(2)=time2-time1
       
      !  open(1, file = 'mps_data.dat')  
      !  do i=1,size(rmlexp)
      !     write(1,*) rmlexp(i)
      !  end do  
      !  close(1)    

cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 3 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
       if(ifprint.ge.1)
     $   call prinf('=== Step 3 (mp to loc+formta+mpeval) ===*',i,0)
c      ... step 3, convert multipole expansions into local expansions
       call cpu_time(time1)
C$     time1=omp_get_wtime()
cc     zero out mexp
       do k=1,6
         do i=1,nboxes
           do j=1,nexptotp
             do idim=1,nd
               mexp(idim,j,i,k) = 0.0d0
             enddo
           enddo
         enddo
       enddo     
c      init uall,dall,...,etc arrays
       allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
       allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
       allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
       allocate(s3478(24,nthd))
       allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
       allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
       allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
       allocate(e1(5,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd))
       allocate(w2(5,nthd),w4(5,nthd),w6(5,nthd),w8(5,nthd))
       allocate(iboxsubcenters(3,8,nthd))
       allocate(iboxfl(2,8,nthd))
c      figure out allocations needed for iboxsrcind
c      and so on
       nmaxt = 0
       do ibox=1,nboxes
         if(nlist3(ibox).gt.0) then
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = iend-istart+1
           if(npts.gt.nmaxt) nmaxt = npts
         endif
       enddo
       allocate(iboxsrcind(nmaxt,nthd))
       allocate(iboxisort(nmaxt,nthd))
       allocate(iboxisort_tmp(nmaxt,nthd))    

       do ilev=2,nlevels
         allocate(iboxlexp(nd*(nterms(ilev)+1)*
     1            (2*nterms(ilev)+1),8,nthd))

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = scales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            istart = isrcse(1,ibox) 
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            if(npts.gt.0) then
c            rescale the multipole expansion
                call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),
     1                 rscpow,tmp(1,0,-nmax,ithd))
cc                process up down for current box
                call mpoletoexp(nd,tmp(1,0,-nmax,ithd),nterms(ilev),
     1              nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)
                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,1),fexpe,fexpo)
                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,2),fexpe,fexpo)
cc                process north-south for current box
                call rotztoy(nd,nterms(ilev),tmp(1,0,-nmax,ithd),
     1              mptmp(1,ithd),rdminus)
                call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev),
     1              nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)
                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,3),fexpe,fexpo)
                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,4),fexpe,fexpo)
cc                process east-west for current box
                call rotztox(nd,nterms(ilev),tmp(1,0,-nmax,ithd),
     1              mptmp(1,ithd),rdplus)
                call mpoletoexp(nd,mptmp(1,ithd),
     1              nterms(ilev),nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),
     3              mexpf2(1,1,ithd),rlsc)
                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,5),fexpe,fexpo)
                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,6),fexpe,fexpo)

              ! if(ibox.eq.12) then
              !   print *, "boxsize = ", boxsize(ilev)
              !   print *, "nboxes = ", nboxes
              !   print *, "itree(....) = ", itree(ipointer(6)+ibox-1)
              !   print *, "nchild = ", nchild
              !   print *, "scales(ilev) = ", scales(ilev)
              !   print *, "nterms(ilev) = ", nterms(ilev)

              !   open(1, file = 'mps_data.dat')  
              !   do i=1,1
              !     do j=1,405
              !       do k=1,1
              !       write(1,*) real(mexpf1(i,j,k))
              !       enddo
              !     enddo
              !   end do  
              !   close(1)    
                
              !   ! open(1, file = 'mps_data.dat')  
              !   ! do i=1,size(rmlexp)
              !   !     write(1,*) rmlexp(i)
              !   ! end do  
              !   ! close(1)    
              ! endif

            endif
         enddo

        !  if (ilev.eq.3) then
        !   open(1, file = 'mps_data.dat')  
        !   do i=1,6
        !     do k=1,577
        !       do j=1,1627
        !           write(1,*) real(mexp(1,j,k,i))
        !         enddo
        !       enddo
        !     enddo
        !     close(1)    
        !   endif

        ! if (ilev.eq.3) then 
        !  open(1, file = 'mps_data.dat')  
        !   do i=1,1
        !     do j=1,405
        !       do k=1,1
        !       write(1,*) real(mexpf1(i,j,k))
        !       enddo
        !     enddo
        !   end do  
        !   close(1)    
        ! endif

        !  if (ilev.eq.3) then 
        !   open(1, file = 'mps_data.dat')  
        !    do i=1,1
        !      do j=1,1627
        !        do k=1,1
        !        write(1,*) (mexpp1(i,j,k))
        !        enddo
        !      enddo
        !    end do  
        !    close(1)    
        !  endif

      !   if (ilev.eq.3) then
      !   open(1, file = 'mps_data.dat')  
      !  do i=1,1
      !   do j=-nmax,nmax
      !     do k=0,nmax
      !      write(1,*) imag(tmp(1,k,j,i))
      !     enddo
      !   enddo
      !  enddo
      !  close(1)
      !  endif
         
         do i=1,nd
          do j=1,nexptot
            mexpf12(i,j,nthd,1) = mexpf1(i,j,nthd)
            mexpf12(i,j,nthd,2) = mexpf2(i,j,nthd)
          enddo
         enddo
            
cc         loop over parent boxes and ship plane wave
c          expansions to the first child of parent boxes. 
c          The codes are now written from a gathering perspective
c          so the first child of the parent is the one
c          recieving all the local expansions
c          coming from all the lists
         rscpow(0) = 1.0d0
         rtmp = scales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
         do ibox = laddr(1,ilev-1),laddr(2,ilev-1)
           ithd = 0
C$         ithd=omp_get_thread_num()
           ithd = ithd + 1
           npts = 0
           nchild = itree(ipointer(4)+ibox-1)
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = npts + iend-istart+1
           if(npts.gt.0.and.nchild.gt.0) then
            

            if(ibox.eq.12) then
              print *, "ilev = ", ilev
              print *, "boxsize = ", boxsize(ilev)
              print *, "nboxes = ", nboxes
              print *, "itree(....) = ", itree(ipointer(6)+ibox-1)
              print *, "nborsi = ",itree(ipointer(7)+mnbors*(ibox-1)+11)
              print *, "nchild = ", nchild
              print *, "scales(ilev) = ", scales(ilev)
              print *, "nterms(ilev) = ", nterms(ilev)
              print *, "cntlist4 = ", cntlist4
              print *, "list4 = ", list4(1,1:4)
              print *, "mnlist4 = ", mnlist4
              print *, "itree(ipointer(5)) = ", itree(ipointer(5))
              print *, "isep = ", isep
              print *, "nn = ", nn
              print *, "lmptot = ", lmptot
              ! print *, "rscpow = ", rscpow
              ! list4(1,2) = 0
              ! list4(1,4) = 0
              
              open(1, file = 'mps_data.dat')  
              do i=1,size(rmlexp)
                  write(1,*) rmlexp(i)
              end do  
              close(1)    

              open(1, file = 'mps_data2.dat')  
              do i=1,1
                do j=1,405
                  do k=1,1
                   write(1,*) real(mexpf2(i,j,k))
                  enddo
                enddo
              end do  
              close(1)    

                open(1, file = 'mps_data3.dat')  
                do i=1,1
                  do j=1,1627
                    do k=1,1
                    write(1,*) real(mexpp1(i,j,k))
                    enddo
                  enddo
                end do  
                close(1)    

                open(1, file = 'mps_data4.dat')  
                do i=1,27037
                  write(1,*) real(fexpback(i))
                end do  
                close(1)    
                
                open(1, file = 'mps_data5.dat') 
                do i=1,16
                  do j=1,1627
                    write(1,*) real(mexppall(1,j,i,1))
                  enddo
                enddo
                close(1)   

                open(1, file = 'mps_data6.dat')
                do i=1,577
                  write(1,*) (nlist4(i))
                enddo
                close(1)  
                
                open(1, file = 'mps_data7.dat')
                do i=1,577
                  write(1,*) (list4ct(i))
                enddo
                close(1)  

                open(1, file = 'mps_data10.dat')
                do i=1,577
                  write(1,*) (list4(1,i))
                enddo
                close(1)  

                open(1, file = 'mps_data9.dat')
                do i=1,577
                  do j=1,3
                    write(1,*) treecenters(j,i)
                  enddo
                enddo

                open(1, file = 'mps_data11.dat')
                do i=1,6
                  do j=1,1627
                    write(1,*) real(pgboxwexp(1,j,1,i))
                  enddo
                enddo

             endif

              call getpwlistallprocessudnsewexp0(ibox,boxsize(ilev),
     1         nboxes,itree(ipointer(6)+ibox-1),
     2         itree(ipointer(7)+mnbors*(ibox-1)),
     3         nchild,itree(ipointer(5)),treecenters,
     4         isep,nd,ilev,scales(ilev),nterms(ilev),iaddrtmp,
     5         rmlexp,lmptottmp,
     6         rlams,whts,nlams,nfourier,nphysical,
     7         nthmax,nexptot,nexptotp,mexp,
     8         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     9         mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     9         mexppall(1,1,1,ithd),rdplus,rdminus,
     9         xshift,yshift,zshift,fexpback,nn,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4) 

              
            !  if(ibox.eq.10) then
            !   print *, "ilev = ", ilev
            !   print *, "boxsize = ", boxsize(ilev)
            !   print *, "nboxes = ", nboxes
            !   print *, "itree(....) = ", itree(ipointer(6)+ibox-1)
            !   print *, "nchild = ", nchild
            !   print *, "scales(ilev) = ", scales(ilev)
            !   print *, "nterms(ilev) = ", nterms(ilev)

            !   open(1, file = 'mps_data.dat')  
            !   do i=1,1
            !     do j=1,405
            !       do k=1,1
            !        write(1,*) imag(mexpf1(i,j,k))
            !       enddo
            !     enddo
            !   end do  
            !   close(1)    
              
            !   ! open(1, file = 'mps_data.dat')  
            !   ! do i=1,size(rmlexp)
            !   !     write(1,*) rmlexp(i)
            !   ! end do  
            !   ! close(1)    
            !  endif
              

           endif
           if(nlist3(ibox).gt.0.and.npts.gt.0) then
             print *, "this is actually processed..."
             do i=1,8
              call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
             enddo
             
             call  getlist3pwallprocessudnsewexp0(ibox,boxsize(ilev),
     1            nboxes,nlist3(ibox),list3(1,ibox),
     2            isep,treecenters,nd,
     3            nterms(ilev),iboxlexp(1,1,ithd),rlams,
     4            whts,nlams,nfourier,nphysical,nthmax,
     5            nexptot,nexptotp,mexp,
     6            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     7            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     8            mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),
     9            xshift,yshift,zshift,fexpback,nn,rlsc,rscpow,
     9            rdplus,rdminus)

            !  if (ibox.eq.148) then
            !   open(1, file = 'mps_data.dat')
            !   do i=1,8
            !     do j=1,nd*(nterms(ilev)+1)*(2*nterms(ilev)+1)
            !       write(1,*) imag(iboxlexp(j,i,ithd))
            !     enddo
            !   enddo
            !   close(1) 
            !  endif
     
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
             if(npts.gt.0) then
               call subdividebox(cmpolesort(1,istart),npts,
     1              treecenters(1,ibox),boxsize(ilev),
     2              iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3              iboxsubcenters(1,1,ithd))
               do i=istart,iend
                 iboxisort_tmp(i-istart+1,ithd) = i
               enddo
               call ireorderf(1,npts,iboxisort_tmp(1,ithd),
     1              iboxisort(1,ithd),iboxsrcind(1,ithd))
               do i=1,8
                 if(iboxfl(1,i,ithd).gt.0) then
                   jstart=iboxfl(1,i,ithd)
                   jend=iboxfl(2,i,ithd)
                   npts0=jend-jstart+1
                   do j = jstart,jend
                     k = iboxisort(j,ithd)
                     call l3dlocloc(nd,scales(ilev),
     1                    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                    nterms(ilev),
     3                    rmpolesort(k),cmpolesort(1,k),
     4                    localsort(impolesort(k)), mtermssort(k),
     5                    dc,lca)
                   enddo
                 endif
               enddo
             endif
           endif
         enddo 
         deallocate(iboxlexp)
        !  if (ilev.eq.2) then
        !  open(1, file = 'mps_data.dat')
        !   do j=1,lmptottmp  
        !     write(1,*) rmlexp(j)
        !   enddo
        !   close(1) 
        ! endif
       enddo

      !  open(1, file = 'mps_data.dat')
      !  do j=1,lmpole  
      !    write(1,*) real(localsort(j))
      !  enddo
      !  close(1) 
      
        ! open(1, file = 'mps_data.dat')
        ! do j=1,lmptottmp  
        !   write(1,*) rmlexp(j)
        ! enddo
        ! close(1) 

      !  open(1, file = 'mps_data.dat')
      !  do j=1,ltree
      !   write(1,*) itree(j)
      !  enddo
       

       deallocate(iboxsrcind,iboxisort,iboxisort_tmp)
       deallocate(iboxsubcenters,iboxfl)
       deallocate(pgboxwexp)
       deallocate(uall,dall,nall,sall,eall,wall)
       deallocate(u1234,d5678,n1256,s3478)
       deallocate(e1357,w2468,n12,n56,s34,s78)
       deallocate(e13,e57,w24,w68)
       deallocate(e1,e3,e5,e7,w2,w4,w6,w8)
       deallocate(tmp,mptmp)
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(3) = timeinfo(3) + time2-time1
cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 4 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
       if(ifprint.ge.1)
     $   call prinf('=== Step 4 (split loc) ===*',i,0)
       call cpu_time(time1)
C$     time1=omp_get_wtime()
       do ilev = 2,nlevels-1
         do ibox = laddr(1,ilev),laddr(2,ilev)

            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
                  if(jbox.gt.0) then
                     call l3dlocloc(nd,scales(ilev),
     1                treecenters(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),scales(ilev+1),treecenters(1,jbox),
     3                rmlexp(iaddr(2,jbox)),nterms(ilev+1),dc,lca)
                  endif
               enddo
            endif
         enddo       
       enddo
       call cpu_time(time2)
C$     time2=omp_get_wtime()
      !  timeinfo(4) = time2-time1
      !  open(1, file = 'mps_data.dat')
      !   do j=1,lmptottmp  
      !     write(1,*) rmlexp(j)
      !   enddo
      !   close(1) 
cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 5 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
       if(ifprint.ge.1)
     $   call prinf('=== step 5 (LOC to CEN) ===*',i,0)
c      ... step 5, shift leaf box loc exp to mpole center
       call cpu_time(time1)
C$     time1=omp_get_wtime()
       do ilev = 0,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
           nchild = itree(ipointer(4)+ibox-1)
           if(nchild.eq.0) then
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend - istart + 1
             do i = istart, iend
              !  print *, "i value: ", i
               call l3dlocloc(nd, scales(ilev),
     1              treecenters(1,ibox), rmlexp(iaddr(2, ibox)),
     2              nterms(ilev), rmpolesort(i), cmpolesort(1,i),
     3              localsort(impolesort(i)), mtermssort(i),
     4              dc,lca)
             enddo
           endif
         enddo
       enddo
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(5) = time2 - time1
      !  open(1, file = 'mps_data.dat')
      !  do j=1,lmpole  
      !    write(1,*) real(localsort(j))
      !  enddo
      !  close(1) 
cccccc
cccccc used to be insdie lfmm3dmain_mps       
cccccc STEP 6 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
       if(ifprint .ge. 1)
     $   call prinf('=== STEP 6 (direct) =====*',i,0)
       call cpu_time(time1)
C$     time1=omp_get_wtime()
cc     directly evaluate potential at sources and targets
c      due to sources in list1
       do ilev=0,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts0 = iend-istart+1
           do iloc = istart,iend
             do i = 1,nlist1(ibox)
               jbox = list1(i,ibox)
               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)
               npts = jend - jstart+1
               do j = jstart,jend
                 d = (cmpolesort(1,j)-cmpolesort(1,iloc))**2
     1             + (cmpolesort(2,j)-cmpolesort(2,iloc))**2
     2             + (cmpolesort(3,j)-cmpolesort(3,iloc))**2
                 d = sqrt(d)
                 if (d .gt. thresh) then
                   call l3dmploc(nd, rmpolesort(j),
     1                  cmpolesort(1,j),
     2                  mpolesort(impolesort(j)), mtermssort(j),
     3                  rmpolesort(iloc), cmpolesort(1,iloc),
     4                  localsort(impolesort(iloc)),
     5                  mtermssort(iloc),dc,lca)
                 endif
               enddo
             enddo
           enddo
         enddo  
       enddo
      !  open(1, file = 'mps_data.dat')
      !  do j=1,lmpole  
      !    write(1,*) imag(localsort(j))
      !  enddo
      !  close(1) 
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(6) = time2-time1
       if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
       d = 0
       do i = 1,6
         d = d + timeinfo(i)
       enddo
       if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
       timeinfo(6) = time2-time1
       if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
       d = 0
       do i = 1,6
         d = d + timeinfo(i)
       enddo
       if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
       if(ier.ne.0) return 
cccccc       
cccccc used to be insdie lfmm3d_mps, after lfmm3dmain_mps
cccccc sor local expansion ccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
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
      !  open(1, file = 'mps_data.dat')
      !  do j=1,nd*ntot  
      !    write(1,*) real(local(j))
      !  enddo
      !  close(1)      
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


c      save data to file in STEP 0       
      ! open(1, file = 'mps_data.dat')  
      !  do i=1,lmpole
      !   write(1,*) imag(mpolesort(i))
      !  enddo
      !  close(1)

      !  open(1, file = 'mps_data.dat')  
      !  do i=1,nmpole
      !   write(1,*) rmpolesort(i)
      !  end do  
      !  close(1)

      !  open(1, file = 'mps_data.dat')  
      !  do i=1,cntlist4
      !   do j=1,8
      !     do k=1,nd*(nterms(nlevels)+1)*(2*nterms(nlevels)+1)
      !       write(1,*) imag(gboxmexp(k,j,i))
      !     enddo
      !   enddo
      !  end do  
      !  close(1)

      !  open(1, file = 'mps_data.dat')  
      !  do i=1,size(rmlexp)
      !     write(1,*) rmlexp(i)
      !  end do  
      !  close(1)    
       
      !  open(1, file = 'mps_data.dat')  
      !  do i=-nmax,nmax
      !   do j=0,nmax
      !     do k=1,nd
      !       write(1,*) imag(tmp(k,j,i,nthd))
      !     enddo
      !   enddo
      !  end do  
      !  close(1)    

      !  open(1, file = 'mps_data.dat')  
      !  do i=1,nexptotp
      !   do j=1,nd
      !       write(1,*) imag(gboxwexp(j,i,1,1,1))
      !   enddo
      !  end do  
      !  close(1)  
       
      !  open(1, file = 'mps_data.dat')  
      !  do k=1,cntlist4
      !  do i=1,nexptotp
      !   do j=1,nd
      !       write(1,*) real(pgboxwexp(j,i,k,6))
      !   enddo
      !  end do 
      !  end do 
      !  close(1)  

      !  open(1, file = 'mps_data.dat')
      !  do j=-nmax,nmax  
      !    do i=0,nmax
      !      do k=0,nmax
      !        write(1,*) rdminus(k,i,j)
      !      enddo
      !    end do  
      !  enddo
      ! close(1)   

      !  output data into a file 
      !  open(1, file = 'mps_data.dat', status = 'new')  
      !  open(1, file = 'mps_data.dat')  
      !  do i=1,size(rmlexp)
      !     write(1,*) rmlexp(i)
      !  end do  
      !  close(1)   
       
      
      !  open(1, file = 'mps_data.dat')  
      !  do i=1,size(rmlexp)
      !     write(1,*) rmlexp(i)
      !  end do  
      !  close(1)   

      !  open(1, file = 'mps_data.dat')
      !  do j=1,nboxes
      !    write(1,*) list4ct(j)
      !  enddo
      !  close(1) 
