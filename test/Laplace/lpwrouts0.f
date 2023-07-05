c----------------------------------------------------------------
c----- to hide a few variables... right now has more though...
      subroutine getpwlistallprocessudnsewexp02(ibox,bs,
     1           nboxes,nnbors,
     2           nbors,
     3           nchild,ichild,centers,
     4           isep,nd,ilev,rscale,nterms,iaddrtmp,
     5           rmlexp,lmptottmp,
     6           rlams,whts,nlams,nfourier,nphysical,
     7           nthmax,nexptot,nexptotp,mexp,mexpupdown,
     9           mexpupphys,mexpdownphys,
     9           mexppall,rdplus,rdminus,
     9           xs,ys,zs,rlsc,rscpow,fexpback,nn)
c-------------------------------------------------------------------
      implicit none
      integer ibox, i, j
      double precision boxsize,bs
      integer nboxes,nnbors,nbors(nnbors)
      integer nchild, ichild(8,nboxes)
      double precision centers(3,nboxes)
      integer isep
      ! additional input for  processudexp0
      integer nd, ilev, nterms, nlams, nthmax, nexptot, nexptotp, nn
      integer nfourier(nlams), nphysical(nlams)
      integer iaddrtmp(2,nboxes)
      integer *8 iaddr(2,nboxes)
      integer lmptottmp
      double precision rscale
      double precision rmlexp(lmptottmp),rlams(nlams),whts(nlams)
      double complex mexp(nd,nexptotp,nboxes,6)
      ! see if possible to pass in only one variable
      double complex mexpupdown(nd,nexptot,2)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexppall(nd,nexptotp,16)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double complex fexpback(nn)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex , allocatable :: pgboxwexp(:,:,:,:)
      integer cntlist4,list4ct(nboxes),nlist4(nboxes),mnlist4
      integer , allocatable :: list4(:,:)
      ! additional input for processnsexp
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      ! additional input for processewexp
      double precision rdminus(0:nterms,0:nterms,-nterms:nterms)


      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

c      init uall,dall,...,etc arrays
      integer uall(200),dall(200)
      integer nall(120),sall(120)
      integer eall(72),wall(72)
      integer u1234(36),d5678(36)
      integer n1256(24),s3478(24)
      integer e1357(16),w2468(16)
      integer n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(5),e3(5),e5(5),e7(5)
      integer w2(5),w4(5),w6(5),w8(5)
     
      integer nthd,ithd

      nthd = 1
      ithd = 1
      iaddr = int(iaddrtmp,8)

c     since list4 is empty... this part needs to be changed...
      cntlist4 = 0
      list4ct = 0
      nlist4 = 0
      mnlist4 = 0
      allocate(pgboxwexp(nd,nexptotp,cntlist4,6))   
      allocate(list4(mnlist4,nboxes))   
     
      call getpwlistall0(ibox,bs,nboxes,nnbors,nbors,
     1           nchild,ichild,centers,isep,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall,nu1234,u1234,
     3           nd5678,d5678,nn1256,n1256,ns3478,s3478,ne1357,e1357,
     4           nw2468,w2468,nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,
     5           e13,ne57,e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,
     6           ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)

      do i=1,nd
        do j=1,nexptot
          mexpup(i,j) = mexpupdown(i,j,1)
          mexpdown(i,j) = mexpupdown(i,j,2)
        enddo
       enddo
  
      call processudexp0(nd,ibox,ilev,nboxes,centers,
     1         ichild,rscale,bs,
     2         nterms,
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nuall,uall,nu1234,u1234,
     5         ndall,dall,nd5678,d5678,
     6         mexpup,mexpdown,mexpupphys,mexpdownphys,
     8         mexppall(1,1,1),mexppall(1,1,2),
     9         mexppall(1,1,3),mexppall(1,1,4),
     9         xs,ys,zs,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)
      
      call processnsexp0(nd,ibox,ilev,nboxes,centers,
     1         ichild,rscale,bs,
     2         nterms,
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nnall,nall,nn1256,n1256,
     5         nn12,n12,nn56,n56,nsall,sall,
     6         ns3478,s3478,ns34,s34,ns78,s78,
     7         mexpup,mexpdown,mexpupphys,mexpdownphys,
     9         mexppall(1,1,1),mexppall(1,1,2),
     9         mexppall(1,1,3),mexppall(1,1,4),
     9         mexppall(1,1,5),mexppall(1,1,6),
     9         mexppall(1,1,7),mexppall(1,1,8),rdplus,
     9         xs,ys,zs,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)

      call processewexp0(nd,ibox,ilev,nboxes,centers,
     1         ichild,rscale,bs,
     2         nterms,
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         neall,eall,ne1357,e1357,
     5         ne13,e13,ne57,e57,ne1,e1,
     6         ne3,e3,ne5,e5,
     7         ne7,e7,nwall,wall,
     8         nw2468,w2468,
     9         nw24,w24,nw68,w68,
     9         nw2,w2,nw4,w4,nw6,w6,
     9         nw8,w8,
     9         mexpup,mexpdown,mexpupphys,mexpdownphys,
     9         mexppall(1,1,1),mexppall(1,1,2),
     9         mexppall(1,1,3),mexppall(1,1,4),
     9         mexppall(1,1,5),mexppall(1,1,6),
     9         mexppall(1,1,7),mexppall(1,1,8),
     9         mexppall(1,1,9),mexppall(1,1,10),
     9         mexppall(1,1,11),mexppall(1,1,12),
     9         mexppall(1,1,13),mexppall(1,1,14),
     9         mexppall(1,1,15),mexppall(1,1,16),rdminus,
     9         xs,ys,zs,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)

       
      
      end subroutine getpwlistallprocessudnsewexp02
      
c----------------------------------------------------------------
c----- to hide a few variables... right now has more though...
      subroutine getpwlistallprocessudnsewexp0(ibox,bs,
     1           nboxes,nnbors,
     2           nbors,
     3           nchild,ichild,centers,
     4           isep,nd,ilev,rscale,nterms,iaddrtmp,
     5           rmlexp,lmptottmp,
     6           rlams,whts,nlams,nfourier,nphysical,
     7           nthmax,nexptot,nexptotp,mexp,
     8           mexpup,mexpdown,
     9           mexpupphys,mexpdownphys,
     9           mexppall,rdplus,rdminus,
     9           xs,ys,zs,fexpback,nn,rlsc,rscpow,
     9           pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)
c-------------------------------------------------------------------
      implicit none
      integer ibox 
      double precision boxsize,bs
      integer nboxes,nnbors,nbors(nnbors)
      integer nchild, ichild(8,nboxes)
      double precision centers(3,nboxes)
      integer isep
      ! additional input for  processudexp0
      integer nd, ilev, nterms, nlams, nthmax, nexptot, nexptotp, nn
      integer nfourier(nlams), nphysical(nlams)
      integer iaddrtmp(2,nboxes)
      integer *8 iaddr(2,nboxes)
      integer lmptottmp
      double precision rscale
      double precision rmlexp(lmptottmp),rlams(nlams),whts(nlams)
      double complex mexp(nd,nexptotp,nboxes,6)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexppall(nd,nexptotp,16)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double complex fexpback(nn)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)
      integer cntlist4,list4ct(nboxes),nlist4(nboxes),mnlist4
      integer list4(mnlist4,nboxes)
      ! additional input for processnsexp
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      ! additional input for processewexp
      double precision rdminus(0:nterms,0:nterms,-nterms:nterms)


      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

c      init uall,dall,...,etc arrays
      integer uall(200),dall(200)
      integer nall(120),sall(120)
      integer eall(72),wall(72)
      integer u1234(36),d5678(36)
      integer n1256(24),s3478(24)
      integer e1357(16),w2468(16)
      integer n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(5),e3(5),e5(5),e7(5)
      integer w2(5),w4(5),w6(5),w8(5)

      integer nthd,ithd

      nthd = 1
      ithd = 1
      iaddr = int(iaddrtmp,8)
     
      call getpwlistall0(ibox,bs,nboxes,nnbors,nbors,
     1           nchild,ichild,centers,isep,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall,nu1234,u1234,
     3           nd5678,d5678,nn1256,n1256,ns3478,s3478,ne1357,e1357,
     4           nw2468,w2468,nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,
     5           e13,ne57,e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,
     6           ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)
  
      call processudexp0(nd,ibox,ilev,nboxes,centers,
     1         ichild,rscale,bs,
     2         nterms,
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nuall,uall,nu1234,u1234,
     5         ndall,dall,nd5678,d5678,
     6         mexpup,mexpdown,mexpupphys,mexpdownphys,
     8         mexppall(1,1,1),mexppall(1,1,2),
     9         mexppall(1,1,3),mexppall(1,1,4),
     9         xs,ys,zs,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)
      
      call processnsexp0(nd,ibox,ilev,nboxes,centers,
     1         ichild,rscale,bs,
     2         nterms,
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nnall,nall,nn1256,n1256,
     5         nn12,n12,nn56,n56,nsall,sall,
     6         ns3478,s3478,ns34,s34,ns78,s78,
     7         mexpup,mexpdown,mexpupphys,mexpdownphys,
     9         mexppall(1,1,1),mexppall(1,1,2),
     9         mexppall(1,1,3),mexppall(1,1,4),
     9         mexppall(1,1,5),mexppall(1,1,6),
     9         mexppall(1,1,7),mexppall(1,1,8),rdplus,
     9         xs,ys,zs,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)

      call processewexp0(nd,ibox,ilev,nboxes,centers,
     1         ichild,rscale,bs,
     2         nterms,
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         neall,eall,ne1357,e1357,
     5         ne13,e13,ne57,e57,ne1,e1,
     6         ne3,e3,ne5,e5,
     7         ne7,e7,nwall,wall,
     8         nw2468,w2468,
     9         nw24,w24,nw68,w68,
     9         nw2,w2,nw4,w4,nw6,w6,
     9         nw8,w8,
     9         mexpup,mexpdown,mexpupphys,mexpdownphys,
     9         mexppall(1,1,1),mexppall(1,1,2),
     9         mexppall(1,1,3),mexppall(1,1,4),
     9         mexppall(1,1,5),mexppall(1,1,6),
     9         mexppall(1,1,7),mexppall(1,1,8),
     9         mexppall(1,1,9),mexppall(1,1,10),
     9         mexppall(1,1,11),mexppall(1,1,12),
     9         mexppall(1,1,13),mexppall(1,1,14),
     9         mexppall(1,1,15),mexppall(1,1,16),rdminus,
     9         xs,ys,zs,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)

       
      
      end subroutine getpwlistallprocessudnsewexp0


c----------------------------------------------------------------
c----- copied from tree_lr_3d.f...
      subroutine getpwlistall0(ibox,bs,nboxes,nnbors,nbors,
     1           nchild,ichild,centers,isep,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall,nu1234,u1234,
     3           nd5678,d5678,nn1256,n1256,ns3478,s3478,ne1357,e1357,
     4           nw2468,w2468,nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,
     5           e13,ne57,e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,
     6           ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)
c-------------------------------------------------------------------
      implicit none
      integer ibox
      double precision boxsize,bs
      integer nboxes,nnbors,nbors(nnbors)
      integer nchild, ichild(8,nboxes)
      double precision centers(3,nboxes)
      integer isep
      integer nuall,ndall,nnall,nsall,neall,nwall,nu1234
      integer nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
      integer uall(200),dall(200)
      integer nall(120),sall(120)
      integer eall(72),wall(72)
      integer u1234(36),d5678(36)
      integer n1256(24),s3478(24)
      integer e1357(16),w2468(16)
      integer n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(5),e3(5),e5(5),e7(5)
      integer w2(5),w4(5),w6(5),w8(5)

      integer jbox,kbox
      integer c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c16
      integer i,j

      nuall = 0
      ndall = 0
      nnall = 0
      nsall = 0
      neall = 0
      nwall = 0
      nu1234 = 0
      nd5678 = 0
      nn1256 = 0
      ns3478 = 0
      ne1357 = 0
      nw2468 = 0
      nn12 = 0
      nn56 = 0
      ns34 = 0
      ns78 = 0
      ne13 = 0
      ne57 = 0
      nw24 = 0
      nw68 = 0
      ne1 = 0
      ne3 = 0
      ne5 = 0
      ne7 = 0
      nw2 = 0
      nw4 = 0
      nw6 = 0
      nw8 = 0
      do i=1,nnbors
         jbox = nbors(i)
         do j=1,8
            kbox = ichild(j,jbox)
            if(kbox.gt.0) then
               c1 = 0
               c2 = 0
               c3 = 0
               c4 = 0
               c5 = 0
               c6 = 0
               c7 = 0
               c8 = 0
               c9 = 0
               c10 = 0
               c11 = 0
               c12 = 0
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c1 = 1
               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c2 = 1
               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c3 = 1
               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c4 = 1
               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c5 = 1
               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c6 = 1
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c7 = 1
               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c8 = 1
               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c9 = 1
               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c10 = 1
               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c11 = 1
               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c12 = 1
               if(c1.eq.1) then
                  nuall = nuall + 1
                  uall(nuall) = kbox
               endif

               if(c2.eq.1) then
                  ndall = ndall + 1
                  dall(ndall) = kbox
               endif

               if(c3.eq.1.and.c1.eq.0.and.c2.eq.0) then
                  nnall = nnall + 1
                  nall(nnall) = kbox
               endif

               if(c4.eq.1.and.c1.eq.0.and.c2.eq.0) then   
                  nsall = nsall + 1
                  sall(nsall) = kbox
               endif

               if(c5.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1             c4.eq.0) then
                  neall = neall + 1
                  eall(neall) = kbox
               endif

               if(c6.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1            c4.eq.0) then
                  nwall = nwall + 1
                  wall(nwall) = kbox
               endif

               c16 = c1 + c2 + c3 + c4 + c5 +c6
               if(c16.eq.0.and.c7.eq.1) then
                  nu1234 = nu1234 + 1
                  u1234(nu1234) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1) then
                  nd5678 = nd5678 + 1
                  d5678(nd5678) = kbox
               endif

               if(c16.eq.0.and.c9.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  nn1256 = nn1256 + 1
                  n1256(nn1256) = kbox
               endif

               if(c16.eq.0.and.c10.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  ns3478 = ns3478 + 1
                  s3478(ns3478) = kbox
               endif

               if(c16.eq.0.and.c11.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  ne1357 = ne1357 + 1
                  e1357(ne1357) = kbox
               endif

               if(c16.eq.0.and.c12.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  nw2468 = nw2468 + 1
                  w2468(nw2468) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c9.eq.1) then
                  nn12 = nn12 + 1
                  n12(nn12) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c9.eq.1) then
                  nn56 = nn56 + 1
                  n56(nn56) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c10.eq.1) then
                  ns34 = ns34 + 1
                  s34(ns34) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c10.eq.1) then
                  ns78 = ns78 + 1
                  s78(ns78) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne13 = ne13 + 1
                  e13(ne13) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne57 = ne57 + 1
                  e57(ne57) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw24 = nw24 + 1
                  w24(nw24) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw68 = nw68 + 1
                  w68(nw68) = kbox
               endif

               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne1 = ne1 + 1
                  e1(ne1) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne3 = ne3 + 1
                  e3(ne3) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne5 = ne5 + 1
                  e5(ne5) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne7 = ne7 + 1
                  e7(ne7) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw2 = nw2 + 1
                  w2(nw2) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw4 = nw4 + 1
                  w4(nw4) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw6 = nw6 + 1
                  w6(nw6) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw8 = nw8 + 1
                  w8(nw8) = kbox
               endif
            endif
         enddo
      enddo
      

      return
      end

c
c
c--------------------------------------------------------------------
      subroutine processudexp0(nd,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           nu1234,u1234,ndall,dall,nd5678,d5678,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpu5678,mexpdall,
     5           mexpd1234,xs,ys,zs,fexpback,rlsc,rscpow,
     6           pgboxwexp,cntlist4,list4,nlist4s,ilist4,mnlist4)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer idim,nd
      integer ibox,ilev,nboxes,nterms,nlams,nthmax
      integer nphysical(nlams),nfourier(nlams)
      integer *8 iaddr(2,nboxes)
      integer ichild(8,nboxes)
      integer nexptot,nexptotp,nmax
      integer nuall,ndall,nu1234,nd5678
      integer uall(*),dall(*),u1234(*),d5678(*)
      double precision rscale,bs
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:)  
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rmlexp(*),centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpuall(nd,nexptotp),mexpdall(nd,nexptotp)
      double complex mexpd1234(nd,nexptotp),mexpu5678(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)
      integer cntlist4,list4(*),nlist4s(*),ilist4(*),mnlist4
      integer nlist4
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)

c      temp variables
      integer jbox,ctr,ii,jj,i,ix,iy,iz,j,kbox
      double precision rtmp,rtmp2
      double complex ztmp,zmul,ztmp2
     
      double precision ctmp(3)

      allocate(tloc(nd,0:nterms,-nterms:nterms))


      do i=1,nexptotp
        do idim=1,nd
          mexpuall(idim,i) = 0
          mexpdall(idim,i) = 0
          mexpu5678(idim,i) = 0
          mexpd1234(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - bs/2.0d0
      ctmp(2) = centers(2,ibox) - bs/2.0d0
      ctmp(3) = centers(3,ibox) - bs/2.0d0
       
      do i=1,nuall
        jbox = uall(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
          do idim=1,nd
            mexpdall(idim,j) = mexpdall(idim,j) + 
     1          mexp(idim,j,jbox,2)*zmul
          enddo
        enddo
      enddo

      do i=1,nu1234
        jbox = u1234(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

        do j=1,nexptotp
          zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
          do idim=1,nd
            mexpd1234(idim,j) = mexpd1234(idim,j) + 
     1          mexp(idim,j,jbox,2)*zmul
          enddo
        enddo
      enddo


      do i=1,ndall
        jbox = dall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

        do j=1,nexptotp
          zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
          do idim=1,nd
            mexpuall(idim,j) = mexpuall(idim,j) + 
     1          mexp(idim,j,jbox,1)*zmul
          enddo
        enddo
      enddo

      do i=1,nd5678
        jbox = d5678(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

        do j=1,nexptotp
          zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
          do idim=1,nd
            mexpu5678(idim,j) = mexpu5678(idim,j) + 
     1         mexp(idim,j,jbox,1)*zmul
          enddo
        enddo
      enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1

      jbox = ichild(1,ibox)

      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)
            mexpdownphys(idim,i) = mexpdall(idim,i) + mexpd1234(idim,i)
          enddo
        enddo
        
        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 2
      jbox = ichild(2,ibox)
      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*xs(1,i)
            mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1          mexpd1234(idim,i))*xs(-1,i)
          enddo
        enddo
 
        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif
  
c      add contributions due to child 3
      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1          mexpd1234(idim,i))*ys(-1,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 4
      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*ztmp
            mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1         mexpd1234(idim,i))*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 5
      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1.0d0/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1           mexpu5678(idim,i))*zs(1,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)*rtmp
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 6
      jbox = ichild(6,ibox)

      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1         mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)


      endif

c      add contributions due to child 7
 
      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1          mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 8
      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)*xs(1,i)
          ztmp2 = xs(-1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1         mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processnsexp0(nd,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nn1256,n1256,nn12,n12,nn56,n56,
     4           nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexpup,
     5           mexpdown,mexpupphys,mexpdownphys,
     6           mexpnall,mexpn3478,mexpn34,mexpn78,mexpsall,
     7           mexps1256,mexps12,mexps56,rdplus,
     8           xs,ys,zs,fexpback,rlsc,rscpow,
     9           pgboxwexp,cntlist4,list4,nlist4s,ilist4,mnlist4)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer nd
      integer ibox,ilev,nboxes,nterms,nlams,nthmax
      integer nphysical(nlams),nfourier(nlams)
      integer *8 iaddr(2,nboxes)
      integer ichild(8,nboxes)
      integer nexptot,nexptotp,nmax
      integer nnall,nsall,nn1256,ns3478,nn12,nn56,ns34,ns78
      integer nall(*),sall(*),n1256(*),s3478(*)
      integer n12(*),n56(*),s34(*),s78(*)
      double precision rscale,bs
      double complex zk2
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:)
      double complex, allocatable :: tloc2(:,:,:)
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      double precision rmlexp(*),centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpnall(nd,nexptotp),mexpsall(nd,nexptotp)
      double complex mexps1256(nd,nexptotp),mexpn3478(nd,nexptotp)
      double complex mexps12(nd,nexptotp),mexps56(nd,nexptotp)
      double complex mexpn34(nd,nexptotp),mexpn78(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)
      integer cntlist4,list4(*),nlist4s(*),ilist4(*),mnlist4
      integer nlist4
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)

c      temp variables
      integer jbox,ctr,ii,jj,i,ix,iy,iz,j,idim,kbox
      double complex ztmp,zmul,ztmp2
      double precision rtmp,rtmp2
    
      double precision ctmp(3)
      allocate(tloc(nd,0:nterms,-nterms:nterms))
      allocate(tloc2(nd,0:nterms,-nterms:nterms))


      do i=1,nexptotp
        do idim=1,nd
          mexpnall(idim,i) = 0
          mexpsall(idim,i) = 0
          mexpn3478(idim,i) = 0
          mexpn34(idim,i) = 0
          mexpn78(idim,i) = 0
          mexps1256(idim,i) = 0
          mexps12(idim,i) = 0
          mexps56(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - bs/2.0d0
      ctmp(2) = centers(2,ibox) - bs/2.0d0
      ctmp(3) = centers(3,ibox) - bs/2.0d0
       
      do i=1,nnall
        jbox = nall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
           zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
           do idim=1,nd
             mexpsall(idim,j) = mexpsall(idim,j) + 
     1           mexp(idim,j,jbox,4)*zmul
           enddo
        enddo

      enddo

      do i=1,nn1256
        jbox = n1256(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
          do idim=1,nd
            mexps1256(idim,j) = mexps1256(idim,j) + 
     1          mexp(idim,j,jbox,4)*zmul
          enddo
        enddo
      enddo

      do i=1,nn12
        jbox = n12(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
          do idim=1,nd
            mexps12(idim,j) = mexps12(idim,j)+mexp(idim,j,jbox,4)*zmul
          enddo
        enddo
      enddo


      do i=1,nn56
        jbox = n56(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
          do idim=1,nd
            mexps56(idim,j) = mexps56(idim,j)+mexp(idim,j,jbox,4)*zmul
          enddo
        enddo
      enddo


      do i=1,nsall
        jbox = sall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
          do idim=1,nd
            mexpnall(idim,j) = mexpnall(idim,j) + 
     1         mexp(idim,j,jbox,3)*zmul
          enddo
        enddo
      enddo

      do i=1,ns3478
        jbox = s3478(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
          do idim=1,nd
            mexpn3478(idim,j) = mexpn3478(idim,j) + 
     1         mexp(idim,j,jbox,3)*zmul
          enddo
        enddo
      enddo

      do i=1,ns34
        jbox = s34(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
          do idim=1,nd
            mexpn34(idim,j) = mexpn34(idim,j)+mexp(idim,j,jbox,3)*zmul
          enddo
        enddo
      enddo

      do i=1,ns78
        jbox = s78(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
          do idim=1,nd
            mexpn78(idim,j) = mexpn78(idim,j)+mexp(idim,j,jbox,3)*zmul
          enddo
        enddo
      enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
      jbox = ichild(1,ibox)

      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)+mexps1256(idim,i)+ 
     1         mexps12(idim,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 2
      jbox = ichild(2,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1         mexps1256(idim,i) + mexps12(idim,i))*ys(-1,i)      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)


      endif
  
c      add contributions due to child 3
      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn34(idim,i)+
     1          mexpn3478(idim,i))*zs(1,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)*rtmp
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 4
      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn34(idim,i)+
     1         mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 5
      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*xs(1,i)
            mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1          mexps1256(idim,i) + mexps56(idim,i))*xs(-1,i)      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 6
      jbox = ichild(6,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ztmp
            mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1          mexps1256(idim,i) + mexps56(idim,i))*ztmp2      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 7
 
      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn78(idim,i)+
     1          mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 8
      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn78(idim,i)+
     1         mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)


        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processewexp0(nd,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,neall,eall,
     3           ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,ne3,e3,ne5,e5,
     4           ne7,e7,nwall,wall,nw2468,w2468,nw24,w24,nw68,w68,
     5           nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6           mexpup,mexpdown,mexpupphys,mexpdownphys,
     7           mexpeall,mexpe2468,mexpe24,mexpe68,mexpe2,mexpe4,
     8           mexpe6,mexpe8,mexpwall,mexpw1357,mexpw13,mexpw57,
     9           mexpw1,mexpw3,mexpw5,mexpw7,rdminus,
     9           xs,ys,zs,fexpback,rlsc,rscpow,
     6           pgboxwexp,cntlist4,list4,nlist4s,ilist4,mnlist4)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer nd
      integer ibox,ilev,nboxes,nterms,nlams,nthmax
      integer nphysical(nlams),nfourier(nlams)
      integer *8 iaddr(2,nboxes)
      integer ichild(8,nboxes)
      integer nexptot,nexptotp,nmax
      integer neall,nwall,ne1357,nw2468,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
      integer eall(*),wall(*),e1357(*),w2468(*)
      integer e13(*),e57(*),w24(*),w68(*)
      integer e1(*),e3(*),e5(*),e7(*),w2(*),w4(*),w6(*),w8(*)
      double precision rscale,bs
      double complex zk2
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:),tloc2(:,:,:)
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rdminus(0:nterms,0:nterms,-nterms:nterms)
      double precision rmlexp(*),centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpeall(nd,nexptotp),mexpwall(nd,nexptotp)
      double complex mexpw1357(nd,nexptotp),mexpe2468(nd,nexptotp)
      double complex mexpw13(nd,nexptotp),mexpw57(nd,nexptotp)
      double complex mexpe24(nd,nexptotp),mexpe68(nd,nexptotp)
      double complex mexpw1(nd,nexptotp),mexpw3(nd,nexptotp)
      double complex mexpw5(nd,nexptotp),mexpw7(nd,nexptotp)
      double complex mexpe2(nd,nexptotp),mexpe4(nd,nexptotp)
      double complex mexpe6(nd,nexptotp),mexpe8(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)
      integer cntlist4,list4(*),nlist4s(*),ilist4(*),mnlist4
      integer nlist4
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)

c      temp variables
      integer jbox,ctr,ii,jj,i,ix,iy,iz,j,l,idim,kbox
      double complex ztmp,zmul,ztmp2
      double precision rtmp,rtmp2
     
      double precision ctmp(3)

      allocate(tloc(nd,0:nterms,-nterms:nterms))
      allocate(tloc2(nd,0:nterms,-nterms:nterms))


      do i=1,nexptotp
        do idim=1,nd
          mexpeall(idim,i) = 0
          mexpwall(idim,i) = 0
          mexpe2468(idim,i) = 0
          mexpe24(idim,i) = 0
          mexpe68(idim,i) = 0
          mexpe2(idim,i) = 0
          mexpe4(idim,i) = 0
          mexpe6(idim,i) = 0
          mexpe8(idim,i) = 0
          mexpw1357(idim,i) = 0
          mexpw13(idim,i) = 0
          mexpw57(idim,i) = 0
          mexpw1(idim,i) = 0
          mexpw3(idim,i) = 0
          mexpw5(idim,i) = 0
          mexpw7(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - bs/2.0d0
      ctmp(2) = centers(2,ibox) - bs/2.0d0
      ctmp(3) = centers(3,ibox) - bs/2.0d0
       
      do i=1,neall
        jbox = eall(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpwall(idim,j) = mexpwall(idim,j) + 
     1         mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,ne1357
        jbox = e1357(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpw1357(idim,j) = mexpw1357(idim,j) + 
     1         mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,ne13
        jbox = e13(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpw13(idim,j) = mexpw13(idim,j)+mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo


      do i=1,ne57
        jbox = e57(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpw57(idim,j) = mexpw57(idim,j)+mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,ne1
        jbox = e1(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpw1(idim,j) = mexpw1(idim,j) + mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo


      do i=1,ne3
        jbox = e3(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
             mexpw3(idim,j) = mexpw3(idim,j) + mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,ne5
        jbox = e5(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
             mexpw5(idim,j) = mexpw5(idim,j) + mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo


      do i=1,ne7
        jbox = e7(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpw7(idim,j) = mexpw7(idim,j) + mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,nwall
        jbox = wall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpeall(idim,j) = mexpeall(idim,j) + 
     1         mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

      do i=1,nw2468
        jbox = w2468(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe2468(idim,j) = mexpe2468(idim,j) + 
     1          mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

      do i=1,nw24
        jbox = w24(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe24(idim,j) = mexpe24(idim,j)+mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

       

      do i=1,nw68
        jbox = w68(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe68(idim,j) = mexpe68(idim,j)+mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

      do i=1,nw2
        jbox = w2(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe2(idim,j) = mexpe2(idim,j) + mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo


      do i=1,nw4
        jbox = w4(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe4(idim,j) = mexpe4(idim,j) + mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

      do i=1,nw6
        jbox = w6(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe6(idim,j) = mexpe6(idim,j) + mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

      do i=1,nw8
        jbox = w8(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpe8(idim,j) = mexpe8(idim,j) + mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
      jbox = ichild(1,ibox)

      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)
            mexpdownphys(idim,i) = mexpwall(idim,i)+mexpw1357(idim,i)+
     1         mexpw13(idim,i)+mexpw1(idim,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 2
      jbox = ichild(2,ibox)

      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1         mexpe24(idim,i)+mexpe2(idim,i))*zs(1,i)      
            mexpdownphys(idim,i) = mexpwall(idim,i)*rtmp
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif
  
c      add contributions due to child 3
      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1         mexpw13(idim,i)+mexpw3(idim,i))*ys(-1,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 4
      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1         mexpe24(idim,i)+mexpe4(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 5
      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*xs(-1,i)
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1             mexpw57(idim,i)+mexpw5(idim,i))*xs(1,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 6
      jbox = ichild(6,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*zs(1,i)
          ztmp2 = xs(1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1           mexpe68(idim,i)+mexpe6(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 7
 
      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ztmp
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1          mexpw57(idim,i)+mexpw7(idim,i))*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 8
      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)*zs(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1         mexpe68(idim,i)+mexpe8(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

      return
      end
