!
!--------------------------------------------------------------------
!
! A fast multi-particle scattering code, based on the code in
! hfmm3d.f of the Flatiron Institute FMM3D library.
!
! Original skeleton code by Manas Rachh, Leslie Greengard, etc.
! FMPS re-write by Mike O'Neil, 2019
!
! The input is assumed to be a collection of multipole expansions,
! and the output is a collection of local expansions at the same
! locations which take into account the potentials from all other
! multipole expansions.
!
! It is assume that all multipole expansions are well-separated, so
! that even those in LIST 1 can be translated.
!
! We use exp(ikr)/r for the Green's function., without the 1/4\pi
! scaling.
! 
!
!--------------------------------------------------------------------
!

subroutine hfmm3d_mps(nd,eps,zk,nsource,source,ifcharge, &
    charge,ifdipole,dipvec, &
    nmpole, cmpole, rmpole, mterms, mpole, impole, lterms, local, &
    ifpgh,pot,grad,hess,ntarg, &
    targ,ifpghtarg,pottarg,gradtarg,hesstarg)
  ! c-----------------------------------------------------------------------
  ! c   INPUT PARAMETERS:
  ! c
  ! c   nd:    in: integer
  ! c             number of densities
  ! c   
  ! c   eps:   in: double precision
  ! c             requested precision
  ! c
  ! c   zk:    in: double complex
  ! c               helmholtz parameter                
  ! c
  ! c   nsource in: integer  
  ! c                number of sources
  ! c
  ! c   source  in: double precision (3,nsource)
  ! c                source(k,j) is the kth component of the jth
  ! c                source locations
  ! c
  ! c   ifcharge  in: integer  
  ! c             charge computation flag
  ! c              ifcharge = 1   =>  include charge contribution
  ! c                                     otherwise do not
  ! c 
  ! c   charge    in: double complex (nd,nsource) 
  ! c              charge strengths
  ! c
  ! c   ifdipole   in: integer
  ! c              dipole computation flag
  ! c              ifdipole = 1   =>  include dipole contribution
  ! c                                     otherwise do not
  ! c
  ! c   dipvec   in: double precision (nd,3,nsource) 
  ! c              dipole orientation vectors
  ! c
  !     nmpole:  in: integer
  !              number of multipole expansion centers
  !
  !     cmpole:  in: double precision (3,nmpole)
  !              multipole expansion centers
  !
  !     rmpole:  in: double precision (nmpole)
  !              scaling factors for each multipole expansion
  !
  !     mterms:  in: integer (nmpole)
  !              order of the multipole expansions, each expansion
  !              can be of a different order
  !
  !     mpole:   in: double complex (nd,0:mterms,-mterms:mterms,nmpole)
  !              coefficients in the multipole expansions
  !
  !     impole:  in: integer (nmpole)
  !              indexing array for mpole, the ith expansion is at
  !              location mpole(impole(i)) and is of order mterms(i)
  !
  ! c   ifpgh   in: integer
  ! c              flag for evaluating potential/gradient at the sources
  ! c              ifpgh = 1, only potential is evaluated
  ! c              ifpgh = 2, potential and gradients are evaluated
  ! c
  ! c
  ! c   ntarg  in: integer  
  ! c                 number of targs 
  ! c
  ! c   targ  in: double precision (3,ntarg)
  ! c               targ(k,j) is the kth component of the jth
  ! c               targ location
  ! c
  ! c   ifpghtarg   in: integer
  ! c              flag for evaluating potential/gradient at the targs
  ! c              ifpghtarg = 1, only potential is evaluated
  ! c              ifpghtarg = 2, potential and gradient are evaluated
  ! c
  ! c
  ! c     OUTPUT parameters:
  ! c
  !     local:   out: double complex ()
  !              local expansions at each center, due to all incoming
  !              multipole expansions (self is ignored)
  !
  ! c   pot:    out: double complex(nd,nsource) 
  ! c               potential at the source locations
  ! c
  ! c   grad:   out: double complex(nd,3,nsource)
  ! c               gradient at the source locations
  ! c
  ! c   hess    out: double complex(nd,6,nsource)
  ! c               hessian at the source locations
  ! c
  ! c   pottarg:    out: double complex(nd,ntarg) 
  ! c               potential at the targ locations
  ! c
  ! c   gradtarg:   out: double complex(nd,3,ntarg)
  ! c               gradient at the targ locations
  ! c
  ! c   hesstarg    out: double complex(nd,6,ntarg)
  ! c                hessian at the target locations
  
  ! c------------------------------------------------------------------
  
  implicit none

  integer nd

  double complex zk
  double precision eps

  integer :: nmpole, mterms, impole(nmpole)
  double precision :: cmpole(3,nmpole), rmpole(nmpole)
  double complex :: mpole(nd,0:mterms,-mterms:mterms,nmpole)
  integer :: lterms
  double complex :: local(nd,*)
  
  integer ifcharge,ifdipole
  integer ifpgh,ifpghtarg

  integer nsource,ntarg

  double precision source(3,nsource),targ(3,ntarg)
  double complex charge(nd,nsource)

  double complex dipvec(nd,3,nsource)

  double complex pot(nd,nsource),grad(nd,3,nsource), &
      pottarg(nd,3,ntarg), &
      gradtarg(nd,3,ntarg),hess(nd,6,*),hesstarg(nd,6,*)

  !
  ! Tree variables
  !
  integer mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
  integer *8 ltree
  integer nlmax
  integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
  integer *8 ipointer(32)
  integer, allocatable :: itree(:)
  double precision, allocatable :: treecenters(:,:),boxsize(:)

  !
  ! temporary sorted arrays
  !
  integer :: lda
  double precision, allocatable :: sourcesort(:,:),targsort(:,:)
  double precision, allocatable :: radsrc(:)
  double complex, allocatable :: chargesort(:,:)
  double complex, allocatable :: dipvecsort(:,:,:)

  double precision, allocatable :: cmpolesort(:,:)
  double precision, allocatable :: rmpolesort(:)
  double complex, allocatable :: mpolesort(:,:,:,:)
  double complex, allocatable :: localsort(:,:,:,:)
  
  double complex, allocatable :: potsort(:,:),gradsort(:,:,:), &
      hesssort(:,:,:)
  double complex, allocatable :: pottargsort(:,:), &
      gradtargsort(:,:,:),hesstargsort(:,:,:)

  !
  !  temporary fmm arrays
  !
  double precision epsfmm
  integer, allocatable :: nterms(:)
  integer *8, allocatable :: iaddr(:,:)
  double precision, allocatable :: scales(:)
  double precision, allocatable :: rmlexp(:)

  integer lmptemp,nmax
  integer *8 lmptot
  double precision, allocatable :: mptemp(:),mptemp2(:)

  !
  !       temporary variables not used in particle code
  !
  double precision expc(3),scjsort(1),radexp
  double complex texpssort(100)
  double precision expcsort(3),radssort(1)
  integer ntj,nexpc,nadd

  !
  !        other temporary variables
  !
  integer i,iert,ifprint,ilev,idim,ier
  double precision time1,time2,omp_get_wtime,second


  !
  !
  ! ifprint is an internal information printing flag.  Suppressed if
  ! ifprint=0.  Prints timing breakdown and other things if ifprint=1.
  !      
  ifprint=1

  !
  ! figure out tree structure
  !
  ! set criterion for box subdivision
  !
  ! if(eps.ge.0.5d-0) then
  !   ndiv = 300
  ! else if(eps.ge.0.5d-1) then
  !   ndiv = 300
  ! else if(eps.ge.0.5d-2) then
  !   ndiv = 300
  ! else if(eps.ge.0.5d-3) then
  !   ndiv = 300
  ! else if(eps.ge.0.5d-6) then
  !   ndiv = 1000
  ! else if(eps.ge.0.5d-9) then
  !   ndiv = 1000
  ! else if(eps.ge.0.5d-12) then
  !   ndiv = 1000
  ! else if(eps.ge.0.5d-15) then
  !   ndiv = 1000
  ! else
  !   ndiv = nsource+ntarg
  ! endif

  ndiv = 1

  if(ifprint.ge.1) print *, "ndiv =",ndiv
  !stop



  !
  ! set tree flags
  !
  isep = 1
  nlmax = 200
  nlevels = 0
  nboxes = 0
  mhung = 0
  ltree = 0

  nexpc = 0
  radexp = 0
  nadd = 0
  ntj = 0

  idivflag = 0

  mnlist1 = 0
  mnlist2 = 0
  mnlist3 = 0
  mnlist4 = 0
  nbmax = 0

  allocate(radsrc(nsource))
  !$omp parallel do default(shared) private(i)
  do i=1,nsource
    radsrc(i) = 0
  enddo
  !$omp end parallel do   



  !
  ! memory management code for constructing level restricted tree
  !
  iert = 0
  call mklraptreemem(iert,source,nsource,radsrc,targ,ntarg, &
      expc,nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax, &
      nlevels,nboxes,mnbors,mnlist1,mnlist2,mnlist3, &
      mnlist4,mhung,ltree)

  if(ifprint.ge.1) print *, ltree/1.0d9


  if(iert.ne.0) then
    print *, "Error in allocating tree memory"
    stop
  endif


  allocate(itree(ltree))
  allocate(boxsize(0:nlevels))
  allocate(treecenters(3,nboxes))

  !
  ! Call tree code
  !
  call mklraptree(source,nsource,radsrc,targ,ntarg,expc, &
      nexpc,radexp,idivflag,ndiv,isep,mhung,mnbors, &
      mnlist1,mnlist2,mnlist3,mnlist4,nlevels, &
      nboxes,treecenters,boxsize,itree,ltree,ipointer)

  !
  !     Allocate sorted source and target arrays      
  !
  
  allocate(sourcesort(3,nsource))
  allocate(targsort(3,ntarg))
  if(ifcharge.eq.1) allocate(chargesort(nd,nsource))
  if(ifdipole.eq.1) then
    allocate(dipvecsort(nd,3,nsource))
  endif

  if(ifpgh.eq.1) then 
    allocate(potsort(nd,nsource),gradsort(nd,3,1),hesssort(nd,6,1))
  else if(ifpgh.eq.2) then
    allocate(potsort(nd,nsource),gradsort(nd,3,nsource), &
        hesssort(nd,6,1))
  else if(ifpgh.eq.3) then
    allocate(potsort(nd,nsource),gradsort(nd,3,nsource), &
        hesssort(nd,6,nsource))
  else
    allocate(potsort(nd,1),gradsort(nd,3,1),hesssort(nd,6,1))
  endif

  if(ifpghtarg.eq.1) then
    allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,1), &
        hesstargsort(nd,6,1))
  else if(ifpghtarg.eq.2) then
    allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg), &
        hesstargsort(nd,6,1))
  else if(ifpghtarg.eq.3) then
    allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg), &
        hesstargsort(nd,6,ntarg))
  else
    allocate(pottargsort(nd,1),gradtargsort(nd,3,1), &
        hesstargsort(nd,6,1))
  endif


  !
  ! scaling factor for multipole and local expansions at all levels
  !
  allocate(scales(0:nlevels),nterms(0:nlevels))
  do ilev = 0,nlevels
    scales(ilev) = boxsize(ilev)
  enddo

  !
  ! initialize potential and gradient at source locations
  !
  if(ifpgh.eq.1) then
    !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
    do i=1,nsource
      do idim=1,nd
        potsort(idim,i) = 0
      enddo
    enddo
    !C$OMP END PARALLEL DO
  endif

  if(ifpgh.eq.2) then
    !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)

    do i=1,nsource
      do idim=1,nd
        potsort(idim,i) = 0
        gradsort(idim,1,i) = 0
        gradsort(idim,2,i) = 0
        gradsort(idim,3,i) = 0
      enddo
    enddo
    !C$OMP END PARALLEL DO
  endif


  if(ifpgh.eq.3) then
    !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
    do i=1,nsource
      do idim=1,nd
        potsort(idim,i) = 0
        gradsort(idim,1,i) = 0
        gradsort(idim,2,i) = 0
        gradsort(idim,3,i) = 0
        hesssort(idim,1,i) = 0
        hesssort(idim,2,i) = 0
        hesssort(idim,3,i) = 0
        hesssort(idim,4,i) = 0
        hesssort(idim,5,i) = 0
        hesssort(idim,6,i) = 0
      enddo
    enddo
    !C$OMP END PARALLEL DO
  endif



  !c
  !cc       initialize potential and gradient  at targ
  !c        locations
  !c
  if(ifpghtarg.eq.1) then
    !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
    do i=1,ntarg
      do idim=1,nd
        pottargsort(idim,i) = 0
      enddo
    enddo
    !C$OMP END PARALLEL DO
  endif

  if(ifpghtarg.eq.2) then
    !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
    do i=1,ntarg
      do idim=1,nd
        pottargsort(idim,i) = 0
        gradtargsort(idim,1,i) = 0
        gradtargsort(idim,2,i) = 0
        gradtargsort(idim,3,i) = 0
      enddo
    enddo
    !C$OMP END PARALLEL DO
  endif

  if(ifpghtarg.eq.3) then
    !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
    do i=1,ntarg
      do idim=1,nd
        pottargsort(idim,i) = 0
        gradtargsort(idim,1,i) = 0
        gradtargsort(idim,2,i) = 0
        gradtargsort(idim,3,i) = 0
        hesstargsort(idim,1,i) = 0
        hesstargsort(idim,2,i) = 0
        hesstargsort(idim,3,i) = 0
        hesstargsort(idim,4,i) = 0
        hesstargsort(idim,5,i) = 0
        hesstargsort(idim,6,i) = 0
      enddo
    enddo
    !C$OMP END PARALLEL DO
  endif


  !
  ! Compute length of expansions at each level      
  !
  nmax = 0
  do i=0,nlevels
    call h3dterms(boxsize(i),zk,eps,nterms(i))
    if(nterms(i).gt.nmax) nmax = nterms(i)
  enddo


  !       
  ! Multipole and local expansions will be held in workspace in
  ! locations pointed to by array iaddr(2,nboxes).
  !
  ! iiaddr is pointer to iaddr array, itself contained in workspace.
  ! imptemp is pointer for single expansion (dimensioned by nmax)
  !
  ! ... allocate iaddr and temporary arrays
  !

  allocate(iaddr(2,nboxes))
  lmptemp = (nmax+1)*(2*nmax+1)*2*nd
  allocate(mptemp(lmptemp),mptemp2(lmptemp))


  !
  ! reorder multipole expansions, their centers, and rscales
  !
  allocate(cmpolesort(3,nmpole))
  allocate(rmpolesort(nmpole))
  allocate(mpolesort(nd,0:mterms,-mterms:mterms,nmpole))
  
  call dreorderf(3, nmpole, cmpole, cmpolesort, itree(ipointer(5)))
  call dreorderf(1, nmpole, rmpole, rmpolesort, itree(ipointer(5)))

  lda = 2*nd*(mterms+1)*(2*mterms+1)
  call dreorderf(lda, nmpole, mpole, mpolesort, itree(ipointer(5)))


  
  !c
  !cc       reorder sources
  !c
  call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))

  if(ifcharge.eq.1) call dreorderf(2*nd,nsource,charge, &
      chargesort, itree(ipointer(5)))

  if(ifdipole.eq.1) then
    call dreorderf(6*nd,nsource,dipvec,dipvecsort, &
        itree(ipointer(5)))
  endif

  !c
  !cc      reorder targs
  !c
  call dreorderf(3,ntarg,targ,targsort,itree(ipointer(6)))

  !
  ! allocate memory need by multipole, local expansions at all
  ! levels
  !
  ! irmlexp is pointer for workspace need by various fmm routines
  !
  call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
  if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9

  allocate(rmlexp(lmptot),stat=iert)
  if(iert.ne.0) then
    print *, "Cannot allocate mpole expansion workspace"
    print *, "lmptot=", lmptot
    stop
  endif


  !
  ! Memory allocation is complete. 
  ! Call main fmm routine
  !
  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  call hfmm3dmain_mps(nd,eps,zk, &
      nsource,sourcesort,&
      ifcharge,chargesort,&
      ifdipole,dipvecsort,&
      nmpole, cmpolesort, rmpolesort, mterms, mpolesort, &
      ntarg,targsort,nexpc,expcsort,radssort, &
      iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp, &
      itree,ltree,ipointer,isep,ndiv,nlevels, &
      nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4, &
      scales,treecenters,itree(ipointer(1)),nterms, &
      ifpgh,potsort,gradsort,hesssort,ifpghtarg,pottargsort, &
      gradtargsort,hesstargsort,ntj,texpssort,scjsort)

  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  if( ifprint .eq. 1 ) call prin2('time in fmm main=*', &
      time2-time1,1)


  !
  ! de-reorder the output local expansions
  !

  

  if(ifpgh.eq.1) then
    call dreorderi(2*nd,nsource,potsort,pot, &
        itree(ipointer(5)))
  endif
  if(ifpgh.eq.2) then 
    call dreorderi(2*nd,nsource,potsort,pot, &
        itree(ipointer(5)))
    call dreorderi(6*nd,nsource,gradsort,grad, &
        itree(ipointer(5)))
  endif

  if(ifpgh.eq.3) then 
    call dreorderi(2*nd,nsource,potsort,pot, &
        itree(ipointer(5)))
    call dreorderi(6*nd,nsource,gradsort,grad, &
        itree(ipointer(5)))
    call dreorderi(12*nd,nsource,hesssort,hess, &
        itree(ipointer(5)))
  endif


  if(ifpghtarg.eq.1) then
    call dreorderi(2*nd,ntarg,pottargsort,pottarg, &
        itree(ipointer(6)))
  endif

  if(ifpghtarg.eq.2) then
    call dreorderi(2*nd,ntarg,pottargsort,pottarg, &
        itree(ipointer(6)))
    call dreorderi(6*nd,ntarg,gradtargsort,gradtarg, &
        itree(ipointer(6)))
  endif

  if(ifpghtarg.eq.3) then
    call dreorderi(2*nd,ntarg,pottargsort,pottarg, &
        itree(ipointer(6)))
    call dreorderi(6*nd,ntarg,gradtargsort,gradtarg, &
        itree(ipointer(6)))
    call dreorderi(12*nd,ntarg,hesstargsort,hesstarg, &
        itree(ipointer(6)))
  endif


  return
end subroutine hfmm3d_mps






subroutine hfmm3dmain_mps(nd,eps,zk, &
    nsource,sourcesort, &
    ifcharge,chargesort, &
    ifdipole,dipvecsort, &
    nmpole, cmpolesort, rmpolesort, mterms, mpolesort, &
    ntarg,targsort,nexpc,expcsort,radssort, &
    iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp, &
    itree,ltree,ipointer,isep,ndiv,nlevels, &
    nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4, &
    rscales,centers,laddr,nterms,ifpgh,pot,grad,hess, &
    ifpghtarg,pottarg,gradtarg,hesstarg, &
    ntj,jsort,scjsort)
  implicit none

  integer nd
  double precision eps
  integer nsource,ntarg, nexpc
  integer ndiv,nlevels

  integer ifcharge,ifdipole
  integer ifpgh,ifpghtarg

  double complex zk,zk2

  double precision sourcesort(3,nsource)

  double complex chargesort(nd,*)
  double complex dipvecsort(nd,3,*)

  ! input multipole stuff
  integer :: nmpole, mterms
  double precision :: cmpolesort(3,nmpole), rmpolesort(nmpole)
  double complex :: mpolesort(nd,0:mterms,-mterms:mterms,nmpole)

  
  double precision targsort(3,ntarg)

  double complex pot(nd,*),grad(nd,3,*),hess(nd,6,*)
  double complex pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

  integer ntj
  double precision expcsort(3,nexpc)
  double complex jsort(nd,0:ntj,-ntj:ntj,nexpc)


  integer lmptemp
  integer *8 iaddr(2,nboxes), lmptot
  double precision rmlexp(lmptot)
  double precision mptemp(lmptemp)
  double precision mptemp2(lmptemp)

  double precision timeinfo(10)
  double precision centers(3,nboxes)

  !c
  !cc      tree variables
  !c
  integer isep
  integer *8 ltree
  integer laddr(2,0:nlevels)
  integer nterms(0:nlevels)
  integer *8 ipointer(32)
  integer itree(ltree)
  integer nboxes
  double precision rscales(0:nlevels)
  double precision boxsize(0:nlevels)
  !c
  !cc      pw stuff
  !c
  integer nuall,ndall,nnall,nsall,neall,nwall
  integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
  integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
  integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

  integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
  integer u1234(36),d5678(36),n1256(24),s3478(24)
  integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
  integer e13(20),e57(20),w24(20),w68(20)
  integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

  integer ntmax, nexpmax, nlams, nmax, nthmax, nphmax
  double precision, allocatable :: carray(:,:), dc(:,:)
  double precision, allocatable :: rdplus(:,:,:)
  double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
  double precision, allocatable :: rdmsq3(:,:,:)
  double complex, allocatable :: rdminus2(:,:,:),zeyep(:)
  double complex, allocatable :: rdplus2(:,:,:)
  double precision, allocatable :: zmone(:)
  integer nn,nnn

  double complex, allocatable :: rlams(:),whts(:)

  double complex, allocatable :: rlsc(:,:,:)
  integer, allocatable :: nfourier(:), nphysical(:)
  integer nexptot, nexptotp
  double complex, allocatable :: xshift(:,:),yshift(:,:),zshift(:,:)

  double complex, allocatable :: fexp(:),fexpback(:)

  double complex, allocatable :: mexp(:,:,:,:)
  double complex, allocatable :: tmp(:,:,:),tmp2(:,:,:)
  double complex, allocatable :: mexpf1(:,:),mexpf2(:,:)
  double complex, allocatable :: mexpp1(:,:),mexpp2(:,:), &
      mexppall(:,:,:)

  double precision, allocatable :: rsc(:)
  double precision r1

  double precision scjsort(nexpc),radssort(nexpc)

  !c     temp variables
  integer i,j,k,l,ii,jj,kk,ll,idim
  integer ibox,jbox,ilev,npts,npts0
  integer nchild,nlist1,nlist2,nlist3,nlist4

  integer istart,iend,istartt,iendt,istarte,iende
  integer istarts,iends
  integer jstart,jend

  integer ifprint,ifwrite

  integer ifhesstarg
  double precision d,time1,time2,omp_get_wtime

  double precision sourcetmp(3)
  double complex chargetmp(nd)

  integer ix,iy,iz
  double precision rtmp
  double complex zmul

  integer nlege, lw7, lused7, itype
  double precision wlege(40000)

  double precision thresh

  integer mnbors,mnlist1, mnlist2,mnlist3,mnlist4
  double complex eye, ztmp,zmult
  double precision alphaj
  integer ctr,ifinit2
  double precision, allocatable :: xnodes(:),wts(:)
  double precision radius
  integer nquad2
  integer maX_nodes
  double precision pi

  integer istart0,istart1,istartm1,nprin
  double precision :: rtmp1,rtmp2,rtmp3,rtmp4, done
  double complex ima

  integer *8 bigint
  integer iert
  data ima/(0.0d0,1.0d0)/

  integer nlfbox,ier, ifstep2


  ntmax = 1000
  allocate(nfourier(ntmax),nphysical(ntmax))
  allocate(rlams(ntmax),whts(ntmax))

  done = 1
  pi = 4.0d0*atan(done)

  nmax = 0
  do i=0,nlevels
    if(nmax.lt.nterms(i)) nmax = nterms(i)
  enddo

  allocate(rsc(0:nmax))

  !
  ! threshold for computing interactions:
  !
  ! interactions will be ignored for all pairs of sources and targets
  ! which satisfy |r| < thresh where r is the disance between them
  !
  thresh = 2.0d0**(-52)*boxsize(0)

  call prini(6,13)
  write(13,*) 'interaction threshold = ', thresh


  allocate(zeyep(-nmax:nmax),zmone(0:2*nmax))

  zeyep(0) = 1
  zmult = -ima
  do i=1,nmax
    zeyep(i) = zeyep(i-1)*zmult
    zeyep(-i) = zeyep(-i+1)/zmult
  enddo


  zmone(0) = 1
  do i=1,2*nmax
    zmone(i) = -zmone(i-1)
  enddo

  !
  ! ifprint is an internal information printing flag. 
  ! Suppressed if ifprint=0.
  ! Prints timing breakdown and other things if ifprint=1.
  ! Prints timing breakdown, list information, and other things if ifprint=2.
  !       
  ifprint=1

  ! ... set the expansion coefficients to zero
  !C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,idim)
  do i=1,nexpc
    do k=-ntj,ntj
      do j = 0,ntj
        do idim=1,nd
          jsort(idim,j,k,i)=0
        enddo
      enddo
    enddo
  enddo
  !C$OMP END PARALLEL DO


  do i=1,10
    timeinfo(i)=0
  enddo

  max_nodes = 10000
  allocate(xnodes(max_nodes))
  allocate(wts(max_nodes))

  !
  ! ... set all multipole and local expansions to zero
  !
  do ilev = 0,nlevels
    !C$OMP PARALLEL DO DEFAULT(SHARED)
    !C$OMP$PRIVATE(ibox)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
      call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
    enddo
    !C$OMP END PARALLEL DO          
  enddo


  !
  ! set scjsort (these are extra output expansion centers)
  !
  do ilev=0,nlevels
    !C$OMP PARALLEL DO DEFAULT(SHARED)
    !C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
    do ibox=laddr(1,ilev),laddr(2,ilev)
      nchild = itree(ipointer(3)+ibox-1)
      if(nchild.gt.0) then
        istart = itree(ipointer(16)+ibox-1)
        iend = itree(ipointer(17)+ibox-1)
        do i=istart,iend
          scjsort(i) = rscales(ilev)
          radssort(i) = min(radssort(i),boxsize(ilev)/32* &
              sqrt(3.0d0))
        enddo
      endif
    enddo
    !C$OMP END PARALLEL DO
  enddo


  !c    initialize legendre function evaluation routines
  nlege = 100
  lw7 = 40000
  call ylgndrfwini(nlege,wlege,lw7,lused7)


  !
  ! ----- Step 1: Shift incoming multipole expansions to the center
  ! of each leaf-node box -----
  !
  if(ifprint .ge. 1) call prinf('=== STEP 1 (shift mp) ====*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  do ilev=2,nlevels

    nquad2 = nterms(ilev)*2.5
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2,xnodes,wts,ifinit2)

    !!!!!!!!!
    ! this radius has a fudge factor in it, debug in future
    !!!!!!!!!
    radius = boxsize(ilev)/2*sqrt(3.0d0)*1.5d0

    if(ifcharge.eq.1.and.ifdipole.eq.0) then

      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,npts,istart,iend,nchild)
      do ibox=laddr(1,ilev),laddr(2,ilev)

        istart = itree(ipointer(10)+ibox-1)
        iend = itree(ipointer(11)+ibox-1)
        npts = iend-istart+1
        
        nchild = itree(ipointer(3)+ibox-1)

        if((npts.gt.0) .and. (nchild.eq.0)) then
          
          call h3dmpmp(nd, zk, rmpolesort(istart), &
              cmpolesort(1,istart), mpolesort(:,:,:,istart), &
              mterms, rscales(ilev), centers(1,ibox), &
              rmlexp(iaddr(1,ibox)), nterms(ilev), &
              radius, xnodes, wts, nquad2)

          !call h3dformmpc(nd,zk,rscales(ilev), &
          !    sourcesort(1,istart),chargesort(1,istart),npts, &
          !    centers(1,ibox),nterms(ilev), &
          !    rmlexp(iaddr(1,ibox)),wlege,nlege)          

        endif


        
      enddo
      !C$OMP END PARALLEL DO            
    endif

  enddo

  call cpu_time(time2)
  !$ time2 = omp_get_wtime()
  timeinfo(1)=time2-time1



  !
  ! ----- Step 2: List 3 interactions, non-adjacent nearfield of all
  ! boxes, do these translations directly right now
  !

  if(ifprint.ge.1) &
      call prinf('=== STEP 2 (form lo) ===*',i,0)
  call cpu_time(time1)
  !$    time1=omp_get_wtime()


  if(ifcharge.eq.1.and.ifdipole.eq.0) then
    do ilev=2,nlevels

      nquad2 = nterms(ilev)*2.5
      nquad2 = max(6,nquad2)
      ifinit2 = 1
      call legewhts(nquad2,xnodes,wts,ifinit2)
      radius = boxsize(ilev)/2*sqrt(3.0d0)

      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
      !C$OMP$SCHEDULE(DYNAMIC)
      do ibox=laddr(1,ilev),laddr(2,ilev)
        nlist4 = itree(ipointer(26)+ibox-1)        

        do i=1,nlist4
          jbox = itree(ipointer(27)+(ibox-1)*mnlist4+i-1)

          !
          ! Form local expansion for all boxes in list3 of the current
          ! box
          !
          istart = itree(ipointer(10)+jbox-1)
          iend = itree(ipointer(11)+jbox-1)
          npts = iend-istart+1
          if(npts.gt.0) then

            if (npts .gt. 1) then
              print *, 'need to add loop for list 3 mp2loc!!!'
              stop
            end if
            
            
            call h3dmploc(nd, zk, rmpolesort(istart),&
                cmpolesort(:,istart), mpolesort(:,:,:,istart), mterms, &
                rscales(ilev), centers(1,ibox), &
                rmlexp(iaddr(2,ibox)), nterms(ilev), &
                radius, xnodes, wts, nquad2)

            
            !call h3dformtac(nd,zk,rscales(ilev), &
            !    sourcesort(1,istart),chargesort(1,istart),npts, &
            !    centers(1,ibox),nterms(ilev), &
            !    rmlexp(iaddr(2,ibox)),wlege,nlege)
          endif
        enddo
      enddo
      !C$OMP END PARALLEL DO
    enddo
  endif


  call cpu_time(time2)
  !$    time2=omp_get_wtime()
  timeinfo(2)=time2-time1





  !
  ! Step 3: Upward pass, multipole-to-multipole merges
  !
  
  if(ifprint .ge. 1) call prinf('=== STEP 3 (merge mp) ====*',i,0)

  call cpu_time(time1)
  !$    time1=omp_get_wtime()

  do ilev=nlevels-1,0,-1
    nquad2 = nterms(ilev)*2.5
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2,xnodes,wts,ifinit2)
    radius = boxsize(ilev)/2*sqrt(3.0d0)

    !C$OMP PARALLEL DO DEFAULT(SHARED)
    !C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      do i=1,8
        jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
        if(jbox.gt.0) then
          istart = itree(ipointer(10)+jbox-1)
          iend = itree(ipointer(11)+jbox-1)
          npts = iend-istart+1

          if(npts.gt.0) then
            call h3dmpmp(nd,zk,rscales(ilev+1), &
                centers(1,jbox),rmlexp(iaddr(1,jbox)), &
                nterms(ilev+1),rscales(ilev),centers(1,ibox), &
                rmlexp(iaddr(1,ibox)),nterms(ilev), &
                radius,xnodes,wts,nquad2)
          endif
        endif
      enddo
    enddo
    !C$OMP END PARALLEL DO          
  enddo

  call cpu_time(time2)
  !$    time2=omp_get_wtime()
  timeinfo(3)=time2-time1




  !
  ! ----- Step 4: Crossward pass, multipole-to-local -----
  ! Note: This is generally the most expensive component of the FMM
  !
  if(ifprint.ge.1) call prinf('=== Step 4 (mp to loc) ===*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  do ilev = 2,nlevels

    ! load the necessary quadrature for plane waves
    zk2 = zk*boxsize(ilev)
    if(real(zk2).le.16*pi .and. imag(zk2).le.12*pi) then
      ier = 0

      ! get new pw quadrature
      call hwts3e(ier,eps,zk2,rlams,whts,nlams)
      call hnumfour(eps,zk2,nlams,nfourier)
      call hnumphys(eps,zk2,nlams,nphysical)

      nphmax = 0
      nthmax = 0
      nexptotp = 0
      nexptot = 0
      nn = 0
      do i=1,nlams
        nexptotp = nexptotp + nphysical(i)
        nexptot = nexptot + 2*nfourier(i)+1
        nn = nn + nfourier(i)*nphysical(i)
        if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
        if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
      enddo

      allocate(fexp(nn),fexpback(nn))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))
      allocate(rlsc(0:nterms(ilev),0:nterms(ilev),nlams))
      allocate(tmp(nd,0:nterms(ilev),-nterms(ilev):nterms(ilev)))
      allocate(tmp2(nd,0:nterms(ilev),-nterms(ilev):nterms(ilev)))

      allocate(mexpf1(nd,nexptot),mexpf2(nd,nexptot), &
          mexpp1(nd,nexptotp))
      allocate(mexpp2(nd,nexptotp),mexppall(nd,nexptotp,16))


      ! NOTE: there can be some additional memory savings here
      bigint = 0
      bigint = nboxes
      bigint = bigint*6
      bigint = bigint*nexptotp*nd

      if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9


      allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
      if(iert.ne.0) then
        print *, "Cannot allocate pw expansion workspace"
        print *, "bigint=", bigint
        stop
      endif


      nn = nterms(ilev)
      allocate(carray(4*nn+1,4*nn+1))
      allocate(dc(0:4*nn,0:4*nn))
      allocate(rdplus(0:nn,0:nn,-nn:nn))
      allocate(rdminus(0:nn,0:nn,-nn:nn))
      allocate(rdsq3(0:nn,0:nn,-nn:nn))
      allocate(rdmsq3(0:nn,0:nn,-nn:nn))

      
      ! generate rotation matrices and carray
      call getpwrotmat(nn,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)

      call hrlscini(rlsc,nlams,rlams,zk2,nterms(ilev))
      call hmkexps(rlams,nlams,nphysical,nexptotp,zk2,xshift, &
          yshift,zshift)

      call hmkfexp(nlams,nfourier,nphysical,fexp,fexpback)

      !
      ! initialize mexp
      !
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(idim,i,j,k)
      do k=1,6
        do i=1,nboxes
          do j=1,nexptotp
            do idim=1,nd
              mexp(idim,j,i,k) = 0.0d0
            enddo
          enddo
        enddo
      enddo
      !C$OMP END PARALLEL DO    

      !
      ! compute powers of scaling parameter for rescaling the
      ! multipole expansions
      !
      r1 = rscales(ilev)
      rsc(0) = 1.0d0
      do i=1,nterms(ilev)
        rsc(i) = rsc(i-1)*r1
      enddo

      !
      ! create multipole to plane wave expansion for all boxes at this
      ! level
      !
      !C$OMP PARALLEL DO DEFAULT (SHARED)
      !C$OMP$PRIVATE(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,tmp2)
      do ibox = laddr(1,ilev),laddr(2,ilev)
        istart = itree(ipointer(10)+ibox-1)
        iend = itree(ipointer(11)+ibox-1)
        npts = iend - istart+1
        if(npts.gt.0) then

          ! rescale multipole expansion
          call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)), &
              rsc,tmp)

          call hmpoletoexp(nd,tmp,nterms(ilev), &
              nlams,nfourier,nexptot,mexpf1,mexpf2,rlsc) 

          call hftophys(nd,mexpf1,nlams,nfourier,nphysical, &
              mexp(1,1,ibox,1),fexp)           

          call hftophys(nd,mexpf2,nlams,nfourier,nphysical, &
              mexp(1,1,ibox,2),fexp)


          ! form mexpnorth, mexpsouth for current box

          ! Rotate mpole for computing mexpnorth and
          ! mexpsouth
          call rotztoy(nd,nterms(ilev),tmp, &
              tmp2,rdminus)

          call hmpoletoexp(nd,tmp2,nterms(ilev),nlams, &
              nfourier,nexptot,mexpf1,mexpf2,rlsc)

          call hftophys(nd,mexpf1,nlams,nfourier, &
              nphysical,mexp(1,1,ibox,3),fexp)           

          call hftophys(nd,mexpf2,nlams,nfourier, &
              nphysical,mexp(1,1,ibox,4),fexp)   


          ! Rotate mpole for computing mexpeast, mexpwest
          call rotztox(nd,nterms(ilev),tmp, &
              tmp2,rdplus)
          call hmpoletoexp(nd,tmp2,nterms(ilev),nlams, &
              nfourier,nexptot,mexpf1,mexpf2,rlsc)

          call hftophys(nd,mexpf1,nlams,nfourier, &
              nphysical,mexp(1,1,ibox,5),fexp)

          call hftophys(nd,mexpf2,nlams,nfourier, &
              nphysical,mexp(1,1,ibox,6),fexp)           

        endif
      enddo
      !C$OMP END PARALLEL DO       

      !
      ! Loop over parent boxes and ship plane wave expansions to the
      ! first child of parent boxes.
      !
      ! The codes are now written from a gathering perspective so the
      ! first child of the parent is the one recieving all the local
      ! expansions coming from all the lists
      !

      !C$OMP PARALLEL DO DEFAULT (SHARED)
      !C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
      !C$OMP$PRIVATE(mexpf1,mexpf2,mexpp1,mexpp2,mexppall)
      !C$OMP$PRIVATE(nuall,uall,ndall,dall,nnall,nall,nsall,sall)
      !C$OMP$PRIVATE(neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678)
      !C$OMP$PRIVATE(nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468)
      !C$OMP$PRIVATE(nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,e57)
      !C$OMP$PRIVATE(nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7)
      !C$OMP$PRIVATE(nw2,w2,nw4,w4,nw6,w6,nw8,w8)
      do ibox = laddr(1,ilev-1),laddr(2,ilev-1)

        npts = 0

        if(ifpghtarg.gt.0) then
          istart = itree(ipointer(12)+ibox-1)
          iend = itree(ipointer(13)+ibox-1)
          npts = npts + iend-istart+1
        endif

        istart = itree(ipointer(14)+ibox-1)
        iend = itree(ipointer(17)+ibox-1)
        npts = npts + iend-istart+1

        nchild = itree(ipointer(3)+ibox-1)

        if(ifpgh.gt.0) then
          istart = itree(ipointer(10)+ibox-1)
          iend = itree(ipointer(11)+ibox-1)
          npts = npts + iend-istart+1
        endif


        if(npts.gt.0.and.nchild.gt.0) then


          call getpwlistall(ibox,boxsize(ilev),nboxes, &
              itree(ipointer(18)+ibox-1),itree(ipointer(19)+ &
              mnbors*(ibox-1)),nchild,itree(ipointer(4)),centers, &
              isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall, &
              neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678, &
              nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468, &
              nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57, &
              e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7, &
              nw2,w2,nw4,w4,nw6,w6,nw8,w8)


          call hprocessudexp(nd,zk2,ibox,ilev,nboxes,centers, &
              itree(ipointer(4)),rscales(ilev),nterms(ilev), &
              iaddr,rmlexp,rlams,whts, &
              nlams,nfourier,nphysical,nthmax,nexptot,nexptotp, &
              mexp, &
              nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678, &
              mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1), &
              mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4), &
              xshift,yshift,zshift,fexpback,rlsc)


          call hprocessnsexp(nd,zk2,ibox,ilev,nboxes,centers,&
              itree(ipointer(4)),rscales(ilev),nterms(ilev),&
              iaddr,rmlexp,rlams,whts,&
              nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
              nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,&
              ns3478,s3478,ns34,s34,ns78,s78,&
              mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),&
              mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),&
              mexppall(1,1,5),mexppall(1,1,6),mexppall(1,1,7),&
              mexppall(1,1,8),rdplus,xshift,yshift,zshift, &
              fexpback,rlsc)

          call hprocessewexp(nd,zk2,ibox,ilev,nboxes,centers,&
              itree(ipointer(4)),rscales(ilev),nterms(ilev),&
              iaddr,rmlexp,rlams,whts,&
              nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
              neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,&
              ne3,e3,ne5,e5,ne7,e7,nwall,wall,&
              nw2468,w2468,nw24,w24,nw68,w68,&
              nw2,w2,nw4,w4,nw6,w6,nw8,w8,&
              mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),&
              mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),&
              mexppall(1,1,5),mexppall(1,1,6),&
              mexppall(1,1,7),mexppall(1,1,8),mexppall(1,1,9),&
              mexppall(1,1,10),mexppall(1,1,11),mexppall(1,1,12),&
              mexppall(1,1,13),mexppall(1,1,14),mexppall(1,1,15),&
              mexppall(1,1,16),rdminus,xshift,yshift,zshift,&
              fexpback,rlsc)
        endif
      enddo
      !C$OMP END PARALLEL DO        

      deallocate(xshift,yshift,zshift,rlsc,tmp,tmp2)
      deallocate(carray,dc,rdplus,rdminus,rdsq3,rdmsq3)

      deallocate(mexpf1,mexpf2,mexpp1,mexpp2,mexppall,mexp)
      deallocate(fexp,fexpback)

    else
      nquad2 = nterms(ilev)*2.2
      nquad2 = max(6,nquad2)
      ifinit2 = 1
      ier = 0

      call legewhts(nquad2,xnodes,wts,ifinit2)

      radius = boxsize(ilev)/2*sqrt(3.0d0)
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,istart,iend,npts,nlist2,i,jbox)
      do ibox = laddr(1,ilev),laddr(2,ilev)

        npts = 0
        if(ifpghtarg.gt.0) then
          istart = itree(ipointer(12)+ibox-1)
          iend = itree(ipointer(13)+ibox-1)
          npts = npts + iend - istart + 1
        endif

        istart = itree(ipointer(14)+ibox-1)
        iend = itree(ipointer(17)+ibox-1)
        npts = npts + iend-istart+1

        if(ifpgh.gt.0) then
          istart = itree(ipointer(10)+ibox-1)
          iend = itree(ipointer(11)+ibox-1)
          npts = npts + iend-istart+1
        endif


        nlist2 = itree(ipointer(22)+ibox-1)
        if(npts.gt.0) then
          do i =1,nlist2
            jbox = itree(ipointer(23)+mnlist2*(ibox-1)+i-1)

            istart = itree(ipointer(10)+jbox-1)
            iend = itree(ipointer(11)+jbox-1)
            npts = iend-istart+1

            if(npts.gt.0) then
              call h3dmploc(nd,zk,rscales(ilev), &
                  centers(1,jbox), &
                  rmlexp(iaddr(1,jbox)),nterms(ilev), &
                  rscales(ilev),centers(1,ibox), &
                  rmlexp(iaddr(2,ibox)),nterms(ilev), &
                  radius,xnodes,wts,nquad2)
            endif
          enddo
        endif
      enddo
      !C$OMP END PARALLEL DO        
    endif
  enddo
  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(4) = time2-time1

  !
  ! ----- Step 5: Downward pass, local-to-local -----
  !
  if(ifprint.ge.1) call prinf('=== Step 5 (split loc) ===*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  do ilev = 2,nlevels-1

    nquad2 = nterms(ilev)*2
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2,xnodes,wts,ifinit2)
    radius = boxsize(ilev+1)/2*sqrt(3.0d0)

    !C$OMP PARALLEL DO DEFAULT(SHARED)
    !C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
    do ibox = laddr(1,ilev),laddr(2,ilev)

      npts = 0

      if(ifpghtarg.gt.0) then
        istart = itree(ipointer(12)+ibox-1)
        iend = itree(ipointer(13)+ibox-1)
        npts = npts + iend-istart+1
      endif

      istart = itree(ipointer(14)+ibox-1)
      iend = itree(ipointer(17)+ibox-1)
      npts = npts + iend-istart+1

      if(ifpgh.gt.0) then
        istart = itree(ipointer(10)+ibox-1)
        iend = itree(ipointer(11)+ibox-1)
        npts = npts + iend-istart+1
      endif

      if(npts.gt.0) then
        do i=1,8
          jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
          if(jbox.gt.0) then
            call h3dlocloc(nd,zk,rscales(ilev), &
                centers(1,ibox),rmlexp(iaddr(2,ibox)), &
                nterms(ilev),rscales(ilev+1),centers(1,jbox), &
                rmlexp(iaddr(2,jbox)),nterms(ilev+1), &
                radius,xnodes,wts,nquad2)
          endif
        enddo
      endif
    enddo
    !C$OMP END PARALLEL DO         
  enddo
  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(5) = time2-time1



  
  if(ifprint.ge.1) call prinf('=== step 6 (mp eval) ===*',i,0)
  call cpu_time(time1)
  !C$        time1=omp_get_wtime()


  !c
  !cc       shift mutlipole expansions to expansion center
  !c        (Note: this part is not relevant for particle codes.
  !c         It is relevant only for QBX codes)

!!! 1000 continue

  nquad2 = max(6,2*ntj)
  ifinit2 = 1
  call legewhts(nquad2,xnodes,wts,ifinit2)
  do ilev=1,nlevels
    !C$OMP PARALLEL DO DEFAULT(SHARED)
    !C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,jbox)
    !C$OMP$SCHEDULE(DYNAMIC)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      nlist3 = itree(ipointer(24)+ibox-1)

      istart = itree(ipointer(16)+ibox-1)
      iend = itree(ipointer(17)+ibox-1)

      do j=istart,iend
        do i=1,nlist3
          jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
          !c
          !cc                  shift multipole expansion directly from box
          !c                   center to expansion center
          call h3dmploc(nd,zk,rscales(ilev+1), &
              centers(1,jbox), &
              rmlexp(iaddr(1,jbox)),nterms(ilev+1), &
              scjsort(j),expcsort(1,j), &
              jsort(1,0,-ntj,j),ntj, &
              radssort(j),xnodes,wts,nquad2)
        enddo
      enddo
    enddo
    !C$OMP END PARALLEL DO  
  enddo

  !c
  !cc       evaluate multipole expansions at source locations
  !c

  do ilev=1,nlevels
    if(ifpgh.eq.1) then         
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
      !C$OMP$SCHEDULE(DYNAMIC)
      do ibox=laddr(1,ilev),laddr(2,ilev)
        nlist3 = itree(ipointer(24)+ibox-1)
        istart = itree(ipointer(10)+ibox-1)
        iend = itree(ipointer(11)+ibox-1)

        npts = iend-istart+1

        do i=1,nlist3
          jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
          call h3dmpevalp(nd,zk,rscales(ilev+1),centers(1,jbox), &
              rmlexp(iaddr(1,jbox)),nterms(ilev+1), &
              sourcesort(1,istart),npts,pot(1,istart),wlege,nlege, &
              thresh)
        enddo
      enddo
      !C$OMP END PARALLEL DO          
    endif

    if(ifpgh.eq.2) then
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
      !C$OMP$SCHEDULE(DYNAMIC)
      do ibox=laddr(1,ilev),laddr(2,ilev)
        nlist3 = itree(ipointer(24)+ibox-1)
        istart = itree(ipointer(10)+ibox-1)
        iend = itree(ipointer(11)+ibox-1)

        npts = iend-istart+1

        do i=1,nlist3
          jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
          call h3dmpevalg(nd,zk,rscales(ilev+1),centers(1,jbox), &
              rmlexp(iaddr(1,jbox)),nterms(ilev+1), &
              sourcesort(1,istart),npts,pot(1,istart), &
              grad(1,1,istart),wlege,nlege,thresh)
        enddo
      enddo
      !C$OMP END PARALLEL DO          
    endif

    if(ifpghtarg.eq.1) then         
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
      !C$OMP$SCHEDULE(DYNAMIC)
      do ibox=laddr(1,ilev),laddr(2,ilev)
        nlist3 = itree(ipointer(24)+ibox-1)
        istart = itree(ipointer(12)+ibox-1)
        iend = itree(ipointer(13)+ibox-1)

        npts = iend-istart+1

        do i=1,nlist3
          jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
          call h3dmpevalp(nd,zk,rscales(ilev+1),centers(1,jbox), &
              rmlexp(iaddr(1,jbox)),nterms(ilev+1), &
              targsort(1,istart),npts,pottarg(1,istart),wlege,nlege, &
              thresh)
        enddo
      enddo
      !C$OMP END PARALLEL DO          
    endif

    if(ifpghtarg.eq.2) then
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
      !C$OMP$SCHEDULE(DYNAMIC)
      do ibox=laddr(1,ilev),laddr(2,ilev)
        nlist3 = itree(ipointer(24)+ibox-1)
        istart = itree(ipointer(12)+ibox-1)
        iend = itree(ipointer(13)+ibox-1)

        npts = iend-istart+1

        do i=1,nlist3
          jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
          call h3dmpevalg(nd,zk,rscales(ilev+1),centers(1,jbox), &
              rmlexp(iaddr(1,jbox)),nterms(ilev+1), &
              targsort(1,istart),npts,pottarg(1,istart), &
              gradtarg(1,1,istart),wlege,nlege,thresh)
        enddo
      enddo
      !C$OMP END PARALLEL DO
    endif
  enddo

  call cpu_time(time2)
  !C$        time2=omp_get_wtime()
  timeinfo(6) = time2-time1





  
  if(ifprint.ge.1) &
      call prinf('=== step 7 (eval lo) ===*',i,0)

  !c     ... step 7, evaluate all local expansions
  !c

  nquad2 = 2*ntj
  nquad2 = max(6,nquad2)
  ifinit2 = 1

  call legewhts(nquad2,xnodes,wts,ifinit2)
  call cpu_time(time1)
  !C$        time1=omp_get_wtime()
  !C

  !c
  !cc       shift local expansion to local epxanion at expansion centers
  !c        (note: this part is not relevant for particle codes.
  !c        it is relevant only for qbx codes)

  do ilev = 0,nlevels
    !C$OMP PARALLEL DO DEFAULT(SHARED)
    !C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
    !C$OMP$SCHEDULE(DYNAMIC)      
    do ibox = laddr(1,ilev),laddr(2,ilev)
      nchild=itree(ipointer(3)+ibox-1)
      if(nchild.eq.0) then 
        istart = itree(ipointer(16)+ibox-1)
        iend = itree(ipointer(17)+ibox-1)
        do i=istart,iend

          call h3dlocloc(nd,zk,rscales(ilev), &
              centers(1,ibox),rmlexp(iaddr(2,ibox)), &
              nterms(ilev),rscales(ilev),expcsort(1,i), &
              jsort(1,0,-ntj,i),ntj,radssort(i),xnodes,wts, &
              nquad2)
        enddo
      endif
    enddo
    !C$OMP END PARALLEL DO
  enddo

  !c
  !cc        evaluate local expansion at source and target
  !c         locations
  !c
  do ilev = 0,nlevels
    if(ifpgh.eq.1) then
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
      !C$OMP$SCHEDULE(DYNAMIC)      
      do ibox = laddr(1,ilev),laddr(2,ilev)
        nchild=itree(ipointer(3)+ibox-1)
        if(nchild.eq.0) then 
          istart = itree(ipointer(10)+ibox-1)
          iend = itree(ipointer(11)+ibox-1)
          npts = iend-istart+1
          call h3dtaevalp(nd,zk,rscales(ilev),centers(1,ibox), &
              rmlexp(iaddr(2,ibox)),nterms(ilev), &
              sourcesort(1,istart), &
              npts,pot(1,istart),wlege,nlege)
        endif
      enddo
      !C$OMP END PARALLEL DO          
    endif

    if(ifpgh.eq.2) then
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
      !C$OMP$SCHEDULE(DYNAMIC)      
      do ibox = laddr(1,ilev),laddr(2,ilev)
        nchild=itree(ipointer(3)+ibox-1)
        if(nchild.eq.0) then 
          istart = itree(ipointer(10)+ibox-1)
          iend = itree(ipointer(11)+ibox-1)
          npts = iend-istart+1
          call h3dtaevalg(nd,zk,rscales(ilev),centers(1,ibox), &
              rmlexp(iaddr(2,ibox)),nterms(ilev), &
              sourcesort(1,istart), &
              npts,pot(1,istart),grad(1,1,istart),wlege,nlege)
        endif
      enddo
      !C$OMP END PARALLEL DO         
    endif

    if(ifpghtarg.eq.1) then
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
      !C$OMP$SCHEDULE(DYNAMIC)      
      do ibox = laddr(1,ilev),laddr(2,ilev)
        nchild=itree(ipointer(3)+ibox-1)
        if(nchild.eq.0) then 
          istart = itree(ipointer(12)+ibox-1)
          iend = itree(ipointer(13)+ibox-1)
          npts = iend-istart+1
          call h3dtaevalp(nd,zk,rscales(ilev),centers(1,ibox), &
              rmlexp(iaddr(2,ibox)),nterms(ilev), &
              targsort(1,istart), &
              npts,pottarg(1,istart),wlege,nlege)
        endif
      enddo
      !C$OMP END PARALLEL DO         
    endif

    if(ifpghtarg.eq.2) then
      !C$OMP PARALLEL DO DEFAULT(SHARED)
      !C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
      !C$OMP$SCHEDULE(DYNAMIC)      
      do ibox = laddr(1,ilev),laddr(2,ilev)
        nchild=itree(ipointer(3)+ibox-1)
        if(nchild.eq.0) then 
          istart = itree(ipointer(12)+ibox-1)
          iend = itree(ipointer(13)+ibox-1)
          npts = iend-istart+1

          call h3dtaevalg(nd,zk,rscales(ilev),centers(1,ibox), &
              rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
              npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
        endif
      enddo
      !C$OMP END PARALLEL DO         
    endif
  enddo


  call cpu_time(time2)
  !C$        time2=omp_get_wtime()
  timeinfo(7) = time2 - time1


  if(ifprint .ge. 1) call prinf('=== STEP 8 (direct) =====*',i,0)

  call cpu_time(time1)
  !C$        time1=omp_get_wtime()

  !c
  !cc       directly form local expansions for list1 sources
  !c        at expansion centers. 
  !c        (note: this part is not relevant for particle codes.
  !c         It is relevant only for qbx codes)


  do ilev=0,nlevels
    !C$OMP PARALLEL DO DEFAULT(SHARED)     
    !C$OMP$PRIVATE(ibox,istarte,iende,nlist1,i,jbox)
    !C$OMP$PRIVATE(jstart,jend)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      istarte = itree(ipointer(16)+ibox-1)
      iende = itree(ipointer(17)+ibox-1)

      nlist1 = itree(ipointer(20)+ibox-1)

      do i =1,nlist1
        jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)


        jstart = itree(ipointer(10)+jbox-1)
        jend = itree(ipointer(11)+jbox-1)

        call hfmm3dexpc_direct(nd,zk,jstart,jend,istarte, &
            iende,sourcesort,ifcharge,chargesort,ifdipole, &
            dipvecsort,expcsort,jsort,scjsort,ntj, &
            wlege,nlege)
      enddo
    enddo
    !C$OMP END PARALLEL DO
  enddo

  !c
  !cc        directly evaluate potential at sources and targets 
  !c         due to sources in list1

  do ilev=0,nlevels
    !c
    !cc           evaluate at the sources
    !c

    if(ifpgh.eq.1) then
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          npts0 = iends-istarts+1
          nlist1 = itree(ipointer(20)+ibox-1)

          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcp_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart),npts,sourcesort(1,istarts), &
                npts0,pot(1,istarts),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO            
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          npts0 = iends-istarts+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectdp_vec_d(nd,zk,sourcesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
                npts0,pot(1,istarts),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          npts0 = iends-istarts+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcdp_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
                npts0,pot(1,istarts),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif
    endif

    if(ifpgh.eq.2) then
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          npts0 = iends-istarts+1
          nlist1 = itree(ipointer(20)+ibox-1)

          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcg_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart),npts,sourcesort(1,istarts), &
                npts0,pot(1,istarts),grad(1,1,istarts),thresh)   
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          npts0 = iends-istarts+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectdg_vec_d(nd,zk,sourcesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
                npts0,pot(1,istarts),grad(1,1,istarts),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          npts0 = iends-istarts+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcdg_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
                npts0,pot(1,istarts),grad(1,1,istarts),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif
    endif

    if(ifpghtarg.eq.1) then
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          npts0 = iendt-istartt+1
          nlist1 = itree(ipointer(20)+ibox-1)

          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcp_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart),npts,targsort(1,istartt), &
                npts0,pottarg(1,istartt),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          npts0 = iendt-istartt+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectdp_vec_d(nd,zk,sourcesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
                npts0,pottarg(1,istartt),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          npts0 = iendt-istartt+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcdp_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
                npts0,pottarg(1,istartt),thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif
    endif

    if(ifpghtarg.eq.2) then
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          npts0 = iendt-istartt+1
          nlist1 = itree(ipointer(20)+ibox-1)

          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcg_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart),npts,targsort(1,istartt), &
                npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
                thresh)   
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          npts0 = iendt-istartt+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectdg_vec_d(nd,zk,sourcesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
                npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
                thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        !C$OMP PARALLEL DO DEFAULT(SHARED)     
        !C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          npts0 = iendt-istartt+1
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1
            call h3ddirectcdg_vec_d(nd,zk,sourcesort(1,jstart), &
                chargesort(1,jstart), &
                dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
                npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
                thresh)          
          enddo
        enddo
        !C$OMP END PARALLEL DO     
      endif
    endif
  enddo

  call cpu_time(time2)
  !C$        time2=omp_get_wtime()
  timeinfo(8) = time2-time1
  if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,8)
  d = 0
  do i = 1,8
    d = d + timeinfo(i)
  enddo

  if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)


  return
end subroutine hfmm3dmain_mps
