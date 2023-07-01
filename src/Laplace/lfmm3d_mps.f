c
       subroutine lfmm3d_mps(nd, eps ,nmpole, cmpole, rmpole, mterms,
     $                   mpole, impole, local, ier)
c
       implicit none

       integer nd,ier
       double precision eps
       integer nmpole, mterms(nmpole), impole(nmpole)
       double precision :: cmpole(3,nmpole), rmpole(nmpole)
       double complex :: mpole(*)
       double complex :: local(*)

c
cc     tree variables
c
       integer idivflag, ndiv, nboxes, nlevels
       integer nlmax, nlmin, iper, ifunif
       integer *8 ipointer(8),ltree
       integer, allocatable :: itree(:)
       integer, allocatable :: isrcse(:,:),isrc(:)
       integer itarg, ntarg
       double precision :: targ(3)
       double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc     temporary sorted arrays
c
       integer :: lmpole, mt, ilen
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
       integer i,iert,ifprint,ilev,ijk,j,l
       double precision time1,time2,omp_get_wtime

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
       ntarg = 0
       targ(1) = 0
       targ(2) = 0
       targ(3) = 0

c
cc     memory management code for contructing level restricted tree
       call pts_tree_mem(cmpole,nmpole,targ,ntarg,idivflag,ndiv,nlmin,
     1      nlmax,iper,ifunif,nlevels,nboxes,ltree)
      
       allocate(itree(ltree))
       allocate(boxsize(0:nlevels))
       allocate(treecenters(3,nboxes))

c       Call tree code
       call pts_tree_build(cmpole,nmpole,targ,ntarg,idivflag,ndiv,
     1      nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,
     2      treecenters,boxsize)
      

       allocate(isrcse(2,nboxes))
       allocate(isrc(nmpole))

       call pts_tree_sort(nmpole,cmpole,itree,ltree,nboxes,nlevels,
     1      ipointer,treecenters,isrc,isrcse)
      
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
C$     time1=omp_get_wtime()
       call lfmm3dmain_mps(nd,eps,
     $      nmpole,cmpolesort,rmpolesort,mtermssort,mpolesort,
     $      impolesort, localsort,
     $      iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $      itree,ltree,ipointer,ndiv,nlevels,
     $      nboxes,iper,boxsize,treecenters,isrcse,
     $      scales,itree(ipointer(1)),nterms,ier)
       if(ier.ne.0) return

       call cpu_time(time2)
C$     time2=omp_get_wtime()
       if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)

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

       return
       end

c       
c---------------------------------------------------------------
c
       subroutine lfmm3dmain_mps(nd,eps,
     $            nmpole,cmpolesort,rmpolesort,mtermssort,mpolesort,
     $            impolesort,localsort,
     $            iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $            itree,ltree,ipointer,ndiv,nlevels,
     $            nboxes,iper,boxsize,centers,isrcse,
     $            rscales,laddr,nterms,ier)
       implicit none

       integer nd,ndiv,nlevels
       integer ier
       double precision eps
       integer nmpole, mtermssort(nmpole)
       double precision :: cmpolesort(3,nmpole), rmpolesort(nmpole)
       double complex :: mpolesort(*),localsort(*)
       integer :: impolesort(nmpole)
       integer nboxes

       integer *8 iaddr(2,nboxes), lmptot
       integer lmptemp
       double precision rmlexp(lmptot)
       double precision mptemp(lmptemp)
       double precision mptemp2(lmptemp)

       double precision thresh

       double precision timeinfo(6)
       double precision centers(3,nboxes)

       integer isep,iper
       integer laddr(2,0:nlevels)
       integer nterms(0:nlevels)
       integer *8 ipointer(8),ltree
       integer itree(ltree)
       double precision rscales(0:nlevels)
       double precision boxsize(0:nlevels)
       integer isrcse(2,nboxes)
       integer, allocatable :: nlist1(:),list1(:,:)
       integer, allocatable :: nlist2(:),list2(:,:)
       integer, allocatable :: nlist3(:),list3(:,:)
       integer, allocatable :: nlist4(:),list4(:,:)

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
       integer i,j,k,l,ii,jj,kk,ll,m,idim,iloc
       integer ibox,jbox,ilev,npts,npts0
       integer nchild

       integer istart,iend
       integer jstart,jend

       integer ifprint

       double precision d,time1,time2,omp_get_wtime

c      PW variables
       integer nlams, nmax, nthmax, nphmax, nmaxt
       integer lca
       double precision, allocatable :: carray(:,:), dc(:,:)
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

       integer nlege, lw7, lused7
       double precision wlege(40000)

       integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
       integer nn
       double precision, allocatable :: rscpow(:)
       double precision pi

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
       integer iert

       integer nthd,ithd
       integer omp_get_max_threads,omp_get_thread_num
       nthd = 1
C$     nthd=omp_get_max_threads()

       pi = 4.0d0*atan(1.0d0)
       thresh = 2.0d0**(-51)*boxsize(0)


c      ifprint is an internal information printing flag.
c      Suppressed if ifprint=0.
c      Prints timing breakdown and other things if ifprint=1.
c      Prints timing breakdown, list information,
c      and other things if ifprint=2.
c       
       ifprint=1
      
c
c      initialize various tree lists
c
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 27

       isep = 1
      
       call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1      centers,itree(ipointer(3)),itree(ipointer(4)),
     2      itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2      itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)

       allocate(list1(mnlist1,nboxes),nlist1(nboxes))
       allocate(list2(mnlist2,nboxes),nlist2(nboxes))
       allocate(list3(mnlist3,nboxes),nlist3(nboxes))
       allocate(list4(mnlist4,nboxes),nlist4(nboxes))

       call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1      centers,itree(ipointer(3)),itree(ipointer(4)),
     2      itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3      itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4      mnlist2,list2,nlist3,mnlist3,list3,
     4      nlist4,mnlist4,list4)
      

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

c      generate rotation matrices and carray
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

       allocate(fexpe(nn),fexpo(nn),fexpback(nn))
       allocate(tmp(nd,0:nmax,-nmax:nmax,nthd))
       allocate(mptmp(lmptemp,nthd))

       allocate(xshift(-5:5,nexptotp))
       allocate(yshift(-5:5,nexptotp))
       allocate(zshift(5,nexptotp))

       allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd),
     1          mexpp1(nd,nexptotp,nthd))
       allocate(mexpp2(nd,nexptotp,nthd),mexppall(nd,nexptotp,16,nthd))

c
cc     NOTE: there can be some memory savings here
c
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

c
c      set all multipole and local expansions to zero
c

       do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
         do ibox=laddr(1,ilev),laddr(2,ilev)
           call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
           call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO        
       enddo


c      initialize legendre function evaluation routines
       nlege = 100
       lw7 = 40000
       call ylgndrfwini(nlege,wlege,lw7,lused7)

c
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
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
C$OMP$REDUCTION(max:nmaxt)
       do ibox=1,nboxes
         if(list4ct(ibox).gt.0) then
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = iend-istart+1
           if(npts.gt.nmaxt) nmaxt = npts
         endif
       enddo
C$OMP END PARALLEL DO

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
         rtmp = rscales(ilev+1)/boxsize(ilev+1)
         do i=1,nterms(ilev+1)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i,j,k)
C$OMP$PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            if(list4ct(ibox).gt.0) then
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1

              if(npts.gt.0) then
                call subdividebox(cmpolesort(1,istart),npts,
     1               centers(1,ibox),boxsize(ilev+1),
     2               gboxind(1,ithd),gboxfl(1,1,ithd),
     3               gboxsubcenters(1,1,ithd))

                do i=istart,iend
                  gboxisort_tmp(i-istart+1,ithd) = i
                enddo

                call ireorderf(1,npts,gboxisort_tmp(1,ithd),
     1               gboxisort(1,ithd),gboxind(1,ithd))

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
     2                     rscales(ilev+1),gboxsubcenters(1,i,ithd),
     3                     gboxmexp(1,i,jbox),nterms(ilev+1),dc,lca)
                    enddo

                    call l3dmpmp(nd,rscales(ilev+1),
     1                   gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox),
     2                   nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3                   rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox),
     1                   rscpow,tmp(1,0,-nmax,ithd))
c
cc                process up down for current box
c
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
c
cc                process north-south for current box
c
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
c
cc                process east-west for current box
c
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
C$OMP END PARALLEL DO
       enddo
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
c

c
c
       if(ifprint .ge. 1)
     $   call prinf('=== STEP 1 (shift mp) ====*',i,0)
         call cpu_time(time1)
C$       time1=omp_get_wtime()
c
c       step 1, shift incoming multipole expansion to the center
c       of each leaf-node box

       do ilev=2,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,npts,istart,iend,nchild,i)
         do ibox=laddr(1,ilev),laddr(2,ilev)
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = iend-istart+1
           nchild = itree(ipointer(4)+ibox-1)

           if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
             do i = istart,iend
               call l3dmpmp(nd,rmpolesort(i),cmpolesort(1,i),
     1              mpolesort(impolesort(i)),mtermssort(i),
     2              rscales(ilev),centers(1,ibox),
     3              rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
             enddo
           endif
         enddo
C$OMP END PARALLEL DO          
       enddo

       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(1)=time2-time1

c       
       if(ifprint .ge. 1)
     $   call prinf('=== STEP 2 (merge mp) ====*',i,0)
       call cpu_time(time1)
C$     time1=omp_get_wtime()
c
       do ilev=nlevels-1,0,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = isrcse(1,jbox)
                  iend = isrcse(2,jbox)
                  npts = iend-istart+1
                  if(npts.gt.0) then
                     call l3dmpmp(nd,rscales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
                  endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
       enddo

       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(2)=time2-time1

       if(ifprint.ge.1)
     $   call prinf('=== Step 3 (mp to loc+formta+mpeval) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

       call cpu_time(time1)
C$     time1=omp_get_wtime()

c
cc     zero out mexp
c 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,idim)
       do k=1,6
         do i=1,nboxes
           do j=1,nexptotp
             do idim=1,nd
               mexp(idim,j,i,k) = 0.0d0
             enddo
           enddo
         enddo
       enddo
C$OMP END PARALLEL DO      

c      init uall,dall,...,etc arrays
       allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
       allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
       allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
       allocate(s3478(24,nthd))
       allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
       allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
       allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
       allocate(e1(20,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd))
       allocate(w2(5,nthd),w4(5,nthd),w6(5,nthd),w8(5,nthd))
       allocate(iboxsubcenters(3,8,nthd))
       allocate(iboxfl(2,8,nthd))
c
c      figure out allocations needed for iboxsrcind
c      and so on
c
       nmaxt = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
C$OMP$REDUCTION(max:nmaxt)
       do ibox=1,nboxes
         if(nlist3(ibox).gt.0) then
           istart = isrcse(1,ibox)
           iend = isrcse(2,ibox)
           npts = iend-istart+1
           if(npts.gt.nmaxt) nmaxt = npts
         endif
       enddo
C$OMP END PARALLEL DO

       allocate(iboxsrcind(nmaxt,nthd))
       allocate(iboxisort(nmaxt,nthd))
       allocate(iboxisort_tmp(nmaxt,nthd))

       do ilev=2,nlevels
         allocate(iboxlexp(nd*(nterms(ilev)+1)*
     1            (2*nterms(ilev)+1),8,nthd))

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts)
C$OMP$PRIVATE(ithd)
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
c
cc                process up down for current box
c
                call mpoletoexp(nd,tmp(1,0,-nmax,ithd),nterms(ilev),
     1              nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)


                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,1),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,2),fexpe,fexpo)


c
cc                process north-south for current box
c
                call rotztoy(nd,nterms(ilev),tmp(1,0,-nmax,ithd),
     1              mptmp(1,ithd),rdminus)
                call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev),
     1              nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,3),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,4),fexpe,fexpo)

c
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
            endif
         enddo
C$OMP END PARALLEL DO         
c
c
cc         loop over parent boxes and ship plane wave
c          expansions to the first child of parent 
c          boxes. 
c          The codes are now written from a gathering perspective
c
c          so the first child of the parent is the one
c          recieving all the local expansions
c          coming from all the lists
c
c          
c
         rscpow(0) = 1.0d0
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
C$OMP$PRIVATE(nuall,ndall,nnall,nsall)
C$OMP$PRIVATE(neall,nwall,nu1234,nd5678)
C$OMP$PRIVATE(nn1256,ns3478,ne1357,nw2468)
C$OMP$PRIVATE(nn12,nn56,ns34,ns78,ne13,ne57)
C$OMP$PRIVATE(nw24,nw68,ne1,ne3,ne5,ne7)
C$OMP$PRIVATE(nw2,nw4,nw6,nw8)
C$OMP$PRIVATE(npts0,jstart,jend,i,j,k)
C$OMP$PRIVATE(ithd)
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

               call getpwlistall(ibox,boxsize(ilev),nboxes,
     1         itree(ipointer(6)+ibox-1),itree(ipointer(7)+
     2         mnbors*(ibox-1)),nchild,itree(ipointer(5)),centers,
     3         isep,nuall,uall(1,ithd),ndall,dall(1,ithd),
     4         nnall,nall(1,ithd),nsall,sall(1,ithd),
     5         neall,eall(1,ithd),nwall,wall(1,ithd),
     6         nu1234,u1234(1,ithd),nd5678,d5678(1,ithd),
     7         nn1256,n1256(1,ithd),ns3478,s3478(1,ithd),
     8         ne1357,e1357(1,ithd),nw2468,w2468(1,ithd),
     9         nn12,n12(1,ithd),nn56,n56(1,ithd),ns34,s34(1,ithd),
     9         ns78,s78(1,ithd),ne13,e13(1,ithd),ne57,e57(1,ithd),
     9         nw24,w24(1,ithd),nw68,w68(1,ithd),ne1,e1(1,ithd),
     9         ne3,e3(1,ithd),ne5,e5(1,ithd),ne7,e7(1,ithd),
     9         nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),
     9         nw8,w8(1,ithd))

               call processudexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nuall,uall(1,ithd),nu1234,u1234(1,ithd),
     5         ndall,dall(1,ithd),nd5678,d5678(1,ithd),
     6         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     7         mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),
     8         mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9         mexppall(1,1,4,ithd),xshift,
     8         yshift,zshift,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,
     9         nlist4,list4,mnlist4)
               
               call processnsexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nnall,nall(1,ithd),nn1256,n1256(1,ithd),
     5         nn12,n12(1,ithd),nn56,n56(1,ithd),nsall,sall(1,ithd),
     6         ns3478,s3478(1,ithd),ns34,s34(1,ithd),ns78,s78(1,ithd),
     7         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     8         mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),
     9         mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9         mexppall(1,1,4,ithd),
     9         mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),
     9         mexppall(1,1,7,ithd),
     9         mexppall(1,1,8,ithd),rdplus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,
     9         nlist4,list4,mnlist4)

               
               call processewexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         neall,eall(1,ithd),ne1357,e1357(1,ithd),
     5         ne13,e13(1,ithd),ne57,e57(1,ithd),ne1,e1(1,ithd),
     6         ne3,e3(1,ithd),ne5,e5(1,ithd),
     7         ne7,e7(1,ithd),nwall,wall(1,ithd),
     8         nw2468,w2468(1,ithd),
     9         nw24,w24(1,ithd),nw68,w68(1,ithd),
     9         nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),
     9         nw8,w8(1,ithd),
     9         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     9         mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),
     9         mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9         mexppall(1,1,4,ithd),
     9         mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),
     9         mexppall(1,1,7,ithd),mexppall(1,1,8,ithd),
     9         mexppall(1,1,9,ithd),
     9         mexppall(1,1,10,ithd),mexppall(1,1,11,ithd),
     9         mexppall(1,1,12,ithd),
     9         mexppall(1,1,13,ithd),mexppall(1,1,14,ithd),
     9         mexppall(1,1,15,ithd),
     9         mexppall(1,1,16,ithd),rdminus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)


           endif

           if(nlist3(ibox).gt.0.and.npts.gt.0) then
             call getlist3pwlistall(ibox,boxsize(ilev),nboxes,
     1            nlist3(ibox),list3(1,ibox),isep,
     2            centers,nuall,uall(1,ithd),ndall,dall(1,ithd),
     3            nnall,nall(1,ithd),
     4            nsall,sall(1,ithd),neall,eall(1,ithd),
     5            nwall,wall(1,ithd))
             do i=1,8
               call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
             enddo

             call processlist3udexplong(nd,ibox,nboxes,centers,
     1            boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2            whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3            nexptotp,mexp,nuall,uall(1,ithd),ndall,dall(1,ithd),
     4            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     5            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     6            mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),
     7            xshift,yshift,zshift,fexpback,rlsc,rscpow)

             call processlist3nsexplong(nd,ibox,nboxes,centers,
     1            boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2            whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3            nexptotp,mexp,nnall,nall(1,ithd),nsall,sall(1,ithd),
     4            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     5            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     6            mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdplus,
     7            xshift,yshift,zshift,fexpback,rlsc,rscpow)

             call processlist3ewexplong(nd,ibox,nboxes,centers,
     1            boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2            whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3            nexptotp,mexp,neall,eall(1,ithd),nwall,wall(1,ithd),
     4            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     5            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     6            mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdminus,
     7            xshift,yshift,zshift,fexpback,rlsc,rscpow)

             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
             if(npts.gt.0) then
               call subdividebox(cmpolesort(1,istart),npts,
     1              centers(1,ibox),boxsize(ilev),
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
                     call l3dlocloc(nd,rscales(ilev),
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
C$OMP END PARALLEL DO        
         deallocate(iboxlexp)
       enddo

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


       if(ifprint.ge.1)
     $   call prinf('=== Step 4 (split loc) ===*',i,0)

       call cpu_time(time1)
C$     time1=omp_get_wtime()
       do ilev = 2,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
                  if(jbox.gt.0) then
                     call l3dlocloc(nd,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),rscales(ilev+1),centers(1,jbox),
     3                rmlexp(iaddr(2,jbox)),nterms(ilev+1),dc,lca)
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
       enddo
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(4) = time2-time1

       if(ifprint.ge.1)
     $   call prinf('=== step 5 (LOC to CEN) ===*',i,0)

c      ... step 5, shift leaf box loc exp to mpole center
c
       call cpu_time(time1)
C$     time1=omp_get_wtime()

       do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
           nchild = itree(ipointer(4)+ibox-1)
           if(nchild.eq.0) then
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend - istart + 1
             do i = istart, iend
               call l3dlocloc(nd, rscales(ilev),
     1              centers(1,ibox), rmlexp(iaddr(2, ibox)),
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


       if(ifprint .ge. 1)
     $   call prinf('=== STEP 6 (direct) =====*',i,0)
       call cpu_time(time1)
C$     time1=omp_get_wtime()

c
cc     directly evaluate potential at sources and targets
c      due to sources in list1

       do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istart,iend,npts0,i,jbox,jstart,jend,npts,iloc)
C$OMP$PRIVATE(j,d)
C$OMP$SCHEDULE(DYNAMIC)
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
C$OMP END PARALLEL DO     
       enddo
       call cpu_time(time2)
C$     time2=omp_get_wtime()
       timeinfo(6) = time2-time1
       if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
       d = 0
       do i = 1,6
         d = d + timeinfo(i)
       enddo

       if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

       return
       end
