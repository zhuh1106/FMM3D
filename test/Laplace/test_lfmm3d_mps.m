% test lfmm3d_mps
%
%
% variables that need special attention on index...
%    boxsize(0:nlevels)
%    nterms(0:nlevels)
%    scales(0:nlevels)
%    laddr(2,0:nlevels)
%    rscpow(0:nmax)
%    dc(0:4*nmax,0:4*nmax)
%    rdplus(0:nmax,0:nmax,-nmax:nmax)
%    rdminus(0:nmax,0:nmax,-nmax:nmax)
%    rdsq3(0:nmax,0:nmax,-nmax:nmax)
%    rdmsq3(0:nmax,0:nmax,-nmax:nmax)
%    rlsc(0:nmax,0:nmax,nlams)
%    tmp(nd,0:nmax,-nmax:nmax,nthd)
%

       clear
       addpath('../../matlab/')
       dc = zeros(51,51);
       dc = getsqrtbinomialcoeffs_mex(50,dc);
       lca = 4*50;
       done = 1;
       pi = 4*atan(done);
% 
% c      initialize printing routine
%        call prini(6,13)

       nd = 1;
      
      
       n1 = 9;
       ns = n1^3;
       nmpole = ns;
      
       disp(['ns = ', num2str(ns), ' ']);
       disp(['nc = ', num2str(nmpole), ' ']);
        
       nt = 19;
      
       source=zeros(3,ns);
       targ=zeros(3,nt);
       cmpole=zeros(3,nmpole);
       charge=zeros(nd,ns);
       dipvec=zeros(nd,3,ns);
       pot=zeros(nd,ns);
       pot2=zeros(nd,ns);
       grad=zeros(nd,3,ns);
       hess=zeros(nd,6,ns);

       pottarg=zeros(nd,nt);
       gradtarg=zeros(nd,3,nt);
       hesstarg=zeros(nd,6,nt);
       eps = 0.5d-9;
      
       disp('==========================================');
       disp('Testing suite for lfmm3d_mps');
       disp(['Requested precision = ', num2str(eps), ' ']);
       

%      generate sources uniformly in the unit cube 
       h = 1.0d0/(n1+1);
       ijk = 0;
       for i = 1:n1
         for j = 1:n1
           for k = 1:n1
             ijk = ijk + 1;
             source(1,ijk) = h*i;
             source(2,ijk) = h*j;
             source(3,ijk) = h*k;
           end
         end
       end

       dnorm = 0;
       for i=1:ns
         for idim=1:nd
           charge(idim,i) = i/ns;
%            charge(idim,i) = hkrand_mex(0);
           dnorm = dnorm + abs(charge(idim,i))^2;

           dipvec(idim,1,i) = hkrand_mex(0);
           dipvec(idim,2,i) = hkrand_mex(0);
           dipvec(idim,3,i) = hkrand_mex(0);
      
           pot(idim,i) = 0;
           grad(idim,1,i) = 0;
           grad(idim,2,i) = 0;
           grad(idim,3,i) = 0;
         end
       end

       dnorm = sqrt(dnorm);
       for i=1:ns
         for idim = 1:nd
           charge(idim,i) = charge(idim,i)/dnorm;
         end
       end   
        
      
       shift = h/1000;
       disp(['shift = ', num2str(shift), ' ']);
       for i = 1:ns
         cmpole(1,i) = source(1,i) + shift;
         cmpole(2,i) = source(2,i);
         cmpole(3,i) = source(3,i);
       end

%      now form a multipole expansion at each center
       mterms=zeros(nmpole,1);
       impole=zeros(nmpole,1);
      
       ntm = 10;
       ntot = 0;
       for i = 1:nmpole
         mterms(i) = ntm;
         ntot = ntot + (mterms(i)+1)*(2*mterms(i)+1);
       end
      
       mpole = zeros(nd*ntot,1);
      
       impole(1) = 1;
       for i = 1:nmpole-1
         ilen = (mterms(i)+1)*(2*mterms(i)+1);
         impole(i+1) = impole(i) + nd*ilen;
       end
      
       
       nlege = 300;
       lw = 2*(nlege+1)^2; % is this dimension right? = lused
       wlege = zeros(lw,1);

       disp(['before ylgndrfwini, lw = ', num2str(lw), ' ']);
       lused = 0;
       [wlege,lused] = ylgndrfwini_mex(nlege, wlege, lw, lused);
       disp(['before ylgndrfwini, lused = ', num2str(lused), ' ']);

%        call zinitialize(nd*ntot*2, mpole)
       ns1 = 1;
       rscale = 1;
       sc = shift;
       if sc < 1 
         rscale = sc;
       end
       disp(['rscale = ', num2str(rscale), ' ']);

       rmpole = zeros(nmpole,1);

       for i = 1:nmpole
         rmpole(i) = rscale;
         mpoletmp = zeros(nd,mterms(i)+1,2*mterms(i)+1);
         mpoletmp =  l3dformmpc_mex(nd, rscale, source(:,i), charge(:,i), ns1, cmpole(:,i), mterms(i), mpoletmp,wlege, nlege);
         mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) = mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) + mpoletmp(:);
%          keyboard
       end
       %%% there is a 1e-32 difference in mpole between matlab mex and fortran...

%      do the direct calculation
       thresh = 1.0d-15;
       ifcharge = 1;
       ifdipole = 0;
       ifpgh = 1;
       ntarg = nt;
       ifpghtarg = 0;
       ier = 0;
       iper = 0;
%        [pot,grad,hess,pottarg,gradtarg,hesstarg] = lfmm3d_mex(nd, eps, ns, source, ifcharge, charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, ntarg, targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier);
%        disp(['via fmm, potential = ', num2str(pot(1,1:10)), ' ']);
       srcinfo.sources = source;
       srcinfo.charges = charge;
       srcinfo.nd = nd;
       U1 = lfmm3d(eps,srcinfo,1);
       pot = U1.pot;
       disp(['via fmm, potential = ', num2str(U1.pot(1:10)), ' ']);
   
       local = zeros(nd*ntot,1);
       
%      now test source to source, charge,with potentials
       disp([' ']);
       disp(['testing multipoles to locals ']);
       disp(['input: multipole expansions ']);
       disp(['output: local expansions ']);

%%%%%% mps call       
%        local = lfmm3d_mps0_mex(nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);
       % call lfmm3dmps0(nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier)
       disp(['ndiv still needs to be optimized ']);
       ndiv = 1;
%
%c     figure out tree structure
%c     set tree flags
       idivflag = 0;
       nlmax = 51;
       nlevels = 0;
       nboxes = 0;
       ltree = 0;
       nlmin = 0;
       iper = 0;
       ifunif = 0;
%c
       iert = 0;
       ntargtree = 0;
       targtree(1) = 0;
       targtree(2) = 0;
       targtree(3) = 0;
%c     memory management code for contructing level restricted tree
       [nlevels,nboxes,ltree] = pts_tree_mem_mex(cmpole,nmpole,targtree,ntargtree,idivflag,ndiv,nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree);
       itree=zeros(ltree,1);
       boxsize=zeros(nlevels+1,1);
       treecenters=zeros(3,nboxes);
       ipointer=zeros(8,1);
%c     Call tree code
       [itree,ipointer,treecenters,boxsize] = pts_tree_build_mex(cmpole,nmpole,targtree,ntargtree,idivflag,ndiv,nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,treecenters,boxsize);
       isrcse=zeros(2,nboxes);
       isrc=zeros(nmpole,1);
       [isrc,isrcse] = pts_tree_sort_mex(nmpole,cmpole,itree,ltree,nboxes,nlevels,ipointer,treecenters,isrc,isrcse);
%c     end of tree build
       disp(['ltree/1e9 = ', num2str(ltree/1.0d9), ' ']);  
       disp(['nlevels= ', num2str(nlevels), ' ']);  
       nterms=zeros(nlevels+1,1);
%c     Compute length of expansions at each level
       nmax = 0;
       for i=0:nlevels
         nterms(i+1) = l3dterms_mex(eps,nterms(i+1));
         if(nterms(i+1)>nmax), nmax = nterms(i+1); end
       end     
%      Multipole and local expansions will be held in workspace in locations pointed to by array iaddr(2,nboxes).
%      iaddr is pointer to iaddr array, itself contained in workspace.
%      imptemp is pointer for single expansion (dimensioned by nmax)
%      ... allocate iaddr and temporary arrays
       iaddr=zeros(2,nboxes);
       lmptemp=(nmax+1)*(2*nmax+1)*2*nd;
       mptemp=zeros(lmptemp,1); mptemp2=zeros(lmptemp,1);       
%c     reorder pole centers
       cmpolesort=zeros(3,nmpole);
       rmpolesort=zeros(nmpole,1);
       impolesort=zeros(nmpole,1);
       mtermssort=zeros(nmpole,1);
       lmpole = 0;
       for i = 1:nmpole
         lmpole = lmpole + (mterms(i)+1)*(2*mterms(i)+1);
       end
       lmpole = nd * lmpole;
       mpolesort = zeros(lmpole,1);
       cmpolesort = dreorderf_mex(3, nmpole, cmpole, cmpolesort, isrc);
       rmpolesort = dreorderf_mex(1, nmpole, rmpole', rmpolesort', isrc)';
       mtermssort = ireorderf_mex(1, nmpole, mterms', mtermssort', isrc)';
       
       impolesort(1) = 1;
       for i=1:nmpole
         mt = mtermssort(i);
         ilen = (mt+1)*(2*mt+1);
         ijk = 1;
         for j = 1:ilen
           for l = 1:nd
             mpolesort(impolesort(i)+ijk-1)=mpole(impole(isrc(i))+ijk-1);
             ijk = ijk + 1;
           end
         end
         if(i<nmpole) 
           impolesort(i+1) = impolesort(i) + nd*ilen;
         end
       end
%
%      allocate memory need by multipole, local expansions at all levels
%      irmlexp is pointer for workspace need by various fmm routines,
%      ~line 467
       itreetmp = reshape(itree(ipointer(1):ipointer(1)+2*(nlevels+1)-1),2,nlevels+1);
       lmptot = 0;
       % inside mpalloc, (laddr(2,i)-laddr(1,i)+1)*2*nn, for real mlexp rmlexp ?
       [iaddr,lmptot] = mpalloc_mex(nd,itreetmp,iaddr,nlevels,lmptot,nterms,nboxes);
       disp(['lmptot/1.0d9 = ', num2str(lmptot/1.0d9), ' ']);  
       rmlexp=zeros(lmptot,1);
%        rmlexp=ones(lmptot,1);
       localsort=zeros(lmpole,1);
%c     Memory allocation is complete.
%c     scaling factor for multipole and local expansions at all levels
       scales=zeros(nlevels+1);
       for ilev = 0:nlevels
         scales(ilev+1) = boxsize(ilev+1);
       end
       laddr = reshape(itree(ipointer(1):(ipointer(1)+2*(nlevels+1)-1)),2,nlevels+1);
%c     Call main fmm_mps routine
       ier = 0;
       nthd = 1;
       pi = 4.0d0*atan(1.0d0);
       thresh = (2.0d0)^(-51)*boxsize(1);
       ifprint=1;
%c     initialize various tree lists
       mnlist1 = 0;
       mnlist2 = 0;
       mnlist3 = 0;
       mnlist4 = 0;
       mnbors = 27;
       isep = 1;
       [mnlist1,mnlist2,mnlist3,mnlist4] = computemnlists_mex(nlevels,nboxes,reshape(itree(ipointer(1):(ipointer(1)+2*(nlevels+1)-1)),2,nlevels+1),boxsize,...
           treecenters,itree(ipointer(3):(ipointer(3)+nboxes-1)),itree(ipointer(4):(ipointer(4)+nboxes-1)),...
           reshape(itree(ipointer(5):(ipointer(5)+8*nboxes-1)),8,nboxes),isep,itree(ipointer(6):(ipointer(6)+nboxes-1)),mnbors,...
           reshape(itree(ipointer(7):(ipointer(7)+mnbors*nboxes-1)),mnbors,nboxes),iper,mnlist1,mnlist2,mnlist3,mnlist4);
%        computemnlists_mex(nlevels,nboxes,laddr,boxsize,...
%            centers,iparent,nchild,...
%            ichild,isep,nnbors,mnbors,...
%            nbors,iper,mnlist1,mnlist2,mnlist3,mnlist4);
       list1=zeros(mnlist1,nboxes); nlist1=zeros(nboxes,1);
       list2=zeros(mnlist2,nboxes); nlist2=zeros(nboxes,1);
       list3=zeros(mnlist3,nboxes); nlist3=zeros(nboxes,1);
       list4=zeros(mnlist4,nboxes); nlist4=zeros(nboxes,1);
       [nlist1,list1,nlist2,list2,nlist3,list3,nlist4,list4]=computelists_mex(nlevels,nboxes,reshape(itree(ipointer(1):(ipointer(1)+2*(nlevels+1)-1)),2,nlevels+1),boxsize,...
           treecenters,itree(ipointer(3):(ipointer(3)+nboxes-1)),itree(ipointer(4):(ipointer(4)+nboxes-1)),...
           reshape(itree(ipointer(5):(ipointer(5)+8*nboxes-1)),8,nboxes),isep,itree(ipointer(6):(ipointer(6)+nboxes-1)),mnbors,...
           reshape(itree(ipointer(7):(ipointer(7)+mnbors*nboxes-1)),mnbors,nboxes),iper,...
           nlist1,mnlist1,list1,nlist2,mnlist2,list2,nlist3,mnlist3,list3,nlist4,mnlist4,list4);
%      ~ line 548
%c     Initialize routines for plane wave mp loc translation
       if(isep == 1) 
         if(eps>=0.5d-2), nlams = 12; end
         if((eps<0.5d-2)&&(eps>=0.5d-3)), nlams = 12; end
         if((eps<0.5d-3)&&(eps>=0.5d-6)), nlams = 20; end
         if((eps<0.5d-6)&&(eps>=0.5d-9)), nlams = 29; end
         if((eps<0.5d-9)), nlams = 37; end
       end
       if(isep == 2)
         if(eps>=0.5d-3), nlams = 9; end
         if((eps<0.5d-3)&&(eps>=0.5d-6)), nlams = 15; end
         if((eps<0.5d-6)&&(eps>=0.5d-9)), nlams = 22; end
         if((eps<0.5d-9)), nlams = 29; end
       end
       rlams=zeros(nlams,1); whts=zeros(nlams,1);
       nphysical=zeros(nlams,1); nfourier=zeros(nlams,1);
       nmax = 0;
       for i=0:nlevels
         if(nmax<nterms(i+1)), nmax = nterms(i+1); end
       end
       rscpow=zeros(nmax+1,1);
       carray=zeros(4*nmax+1,4*nmax+1);
       dc=zeros(4*nmax+1,4*nmax+1);
       rdplus=zeros(nmax+1,nmax+1,2*nmax+1);
       rdminus=zeros(nmax+1,nmax+1,2*nmax+1);
       rdsq3=zeros(nmax+1,nmax+1,2*nmax+1);
       rdmsq3=zeros(nmax+1,nmax+1,2*nmax+1);
       rlsc=zeros(nmax+1,nmax+1,nlams);
%c     generate rotation matrices and carray
       [carray,rdplus,rdminus,rdsq3,rdmsq3,dc] = getpwrotmat_mex(nmax,carray,rdplus,rdminus,rdsq3,rdmsq3,dc);       
%c     generate rlams and weights (these are the nodes
%c     and weights for the lambda integral)
       [rlams,whts] = vwts_mex(rlams,whts,nlams);
%c     generate the number of fourier modes required to represent the
%c     moment function in fourier space
       nfourier = numthetahalf_mex(nfourier,nlams);
%c     generate the number of fourier modes in physical space
%c     required for the exponential representation 
       nphysical = numthetafour_mex(nphysical,nlams);
%c     Generate powers of lambda for the exponential basis
       rlsc = rlscini_mex(rlsc,nlams,rlams,nmax);
%c     Compute total number of plane waves
       nexptotp = 0;
       nexptot = 0;
       nthmax = 0;
       nphmax = 0;
       nn = 0;
       for i=1:nlams
         nexptot = nexptot + nfourier(i);
         nexptotp = nexptotp + nphysical(i);
         if(nfourier(i)>nthmax), nthmax = nfourier(i); end
         if(nphysical(i)>nphmax), nphmax = nphysical(i); end
         nn = nn + nphysical(i)*nfourier(i);
       end
       fexpe=zeros(nn,1); fexpo=zeros(nn,1); fexpback=zeros(nn,1);
       %%% come back and change nthd to 1...
       tmp=zeros(nd,nmax+1,2*nmax+1,nthd);
       mptmp=zeros(lmptemp,nthd);
       xshift=zeros(11,nexptotp);
       yshift=zeros(11,nexptotp);
       zshift=zeros(5,nexptotp);
       mexpf1=zeros(nd,nexptot,nthd); mexpf2=zeros(nd,nexptot,nthd); mexpp1=zeros(nd,nexptotp,nthd);
       mexpp2=zeros(nd,nexptotp,nthd); mexppall=zeros(nd,nexptotp,16,nthd);
%cc    NOTE: there can be some memory savings here
       bigint = 0;
       bigint = nboxes;
       bigint = bigint*6;
       bigint = bigint*nexptotp*nd;
       list4ct=zeros(nboxes,1);
       ilist4=zeros(nboxes,1);
       for i=1:nboxes
         list4ct(i)=0;
         ilist4(i)=0;
       end
       cntlist4=0;
%c     Precompute table for shifting exponential coefficients in 
%c     physical domain  
       [xshift,yshift,zshift] = mkexps_mex(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift);
%c     Precompute table of exponentials for mapping from
%c     fourier to physical domain
       [fexpe,fexpo,fexpback] = mkfexp_mex(nlams,nfourier,nphysical,fexpe,fexpo,fexpback);
%c     set all multipole and local expansions to zero
       for ilev = 0:nlevels
         for ibox=laddr(1,ilev+1):laddr(2,ilev+1) % laddr(2,0:nlevels)
           mpoleiboxtmp = zeros(nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1),1);
           mpoleiboxtmp = mpzero_mex(nd,mpoleiboxtmp,nterms(ilev+1));
           mpoleiboxtmp = reshape([real(mpoleiboxtmp(:)),imag(mpoleiboxtmp(:))]',[],1);
           tmpidx = 1:2*(nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1));
           rmlexp(iaddr(1,ibox)-1+tmpidx) = mpoleiboxtmp;
           rmlexp(iaddr(2,ibox)-1+tmpidx) = mpoleiboxtmp;
         end      
       end
%c     initialize legendre function evaluation routines
       nlege = 100;
       lw7 = 5*(nlege+1)^2;
       lused7 = 0;
       wlege7 = zeros(lw7,1);
       [wlege7,lused7] = ylgndrfwini_mex(nlege, wlege7, lw7, lused7);
       wlege = wlege7;
%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 0 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc
%c      count number of boxes are in list4
       lca = 4*nmax;
       disp(['=== STEP 0 list4===*']);  
       for ilev=1:nlevels-1
         for ibox=laddr(1,ilev+1):laddr(2,ilev+1)
           if(nlist3(ibox)>0) 
             cntlist4=cntlist4+1;
             list4ct(ibox)=cntlist4;
             ilist4(cntlist4)=ibox;
           end
         end
       end    
       disp(['nboxes: ', num2str(nboxes), ', cntlist4: ', num2str(cntlist4)]);  
       pgboxwexp=zeros(nd,nexptotp,cntlist4,6);
       gboxmexp=zeros(nd*(nterms(nlevels+1)+1)*(2*nterms(nlevels+1)+1),8,cntlist4);
       gboxsubcenters=zeros(3,8,nthd);
       gboxfl=zeros(2,8,nthd);
       nmaxt = 0;
       for ibox=1:nboxes
         if(list4ct(ibox)>0) % this test case never computes below
           istart = isrcse(1,ibox);
           iend = isrcse(2,ibox);
           npts = iend-istart+1;
           if(npts>nmaxt)
             nmaxt = npts;
           end
         end
       end
%        list4ct_f=importdata('mps_data.dat');
       gboxind=zeros(nmaxt,nthd);
       gboxisort=zeros(nmaxt,nthd);
       gboxisort_tmp=zeros(nmaxt,nthd);
       gboxwexp=zeros(nd,nexptotp,6,8,nthd);
%c     note gboxmexp is an array not scalar
%      pgboxwexp=0d0; gboxmexp=0d0;
%c     form mexp for all list4 type box at first ghost box center
       for ilev=1:nlevels-1
         rscpow(1) = 1.0d0/boxsize(ilev+2); % only used inside if(list4ct(ibox)>0)
         rtmp = scales(ilev+2)/boxsize(ilev+2);
         for i=1:nterms(ilev+2) 
           rscpow(i+1) = rscpow(i)*rtmp; 
         end
         for ibox=laddr(1,ilev+1):laddr(2,ilev+1)
            ithd = 0; 
            ithd = ithd + 1;
            if(list4ct(ibox)>0)
%               disp(['list4ct is accessed']);
%               disp(['this is ilev = ' num2str(ilev) ' ']);
%               disp(['this is ibox = ' num2str(ibox) ' ']);
              istart=isrcse(1,ibox); 
              iend=isrcse(2,ibox); 
              npts = iend-istart+1;
              % line 709 onwards
              if(npts>0)
                gboxindi = gboxind(:,ithd);
                gboxfli = gboxfl(:,:,ithd);
                gboxsubcentersi = gboxsubcenters(:,:,ithd);
                [gboxind(:,ithd),gboxfl(:,:,ithd),gboxsubcenters(:,:,ithd)] = ...
                    subdividebox_mex(cmpolesort(:,istart),npts,treecenters(:,ibox),boxsize(ilev+2),gboxindi,gboxfli,gboxsubcentersi);
                for i=istart:iend
                  gboxisort_tmp(i-istart+1,ithd) = i;
                end
                gboxisort_tmpi = gboxisort_tmp(:,ithd);
                gboxisorti = gboxisort(:,ithd);
                gboxindi = gboxind(:,ithd);
                gboxisort(:,ithd) = ireorderf_mex(1, npts, gboxisort_tmpi', gboxisorti', gboxindi)';
                for i=1:8
                  if (gboxfl(1,i,ithd)>0)
                    jstart=gboxfl(1,i,ithd);
                    jend=gboxfl(2,i,ithd);
                    npts0=jend-jstart+1;
                    jbox=list4ct(ibox);
                    for j = jstart:jend % seems like there is always 1 ghost box, why is this the case?
                      k = gboxisort(j,ithd);
                      gboxmexpi = zeros(nd,nterms(ilev+2)+1,2*nterms(ilev+2)+1);
                      gboxmexpi = l3dmpmp_mex(nd,rmpolesort(k),cmpolesort(:,k),...
                          mpolesort(impolesort(k):(impolesort(k)+nd*(mtermssort(k)+1)*(2*mtermssort(k)+1)-1)),mtermssort(k),...
                          scales(ilev+2),gboxsubcenters(:,i,ithd),...
                          gboxmexpi,nterms(ilev+2),dc,lca);
                      gboxmexp(:,i,jbox) = gboxmexp(:,i,jbox) + gboxmexpi(:);
                    end
                    % fortran line 744, from gboxmexp to rmlexp(real).... 
                    cmlexpi = zeros(nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
                    cmlexpi = l3dmpmp_mex(nd,scales(ilev+2),...
                        gboxsubcenters(:,i,ithd),gboxmexp(:,i,jbox),...
                       nterms(ilev+2),scales(ilev+1),treecenters(:,ibox),...
                        cmlexpi,nterms(ilev+1),dc,lca);
                    %%% here pay attention to the variable type, is it real or complex...
                    tmpidx = 1:2*(nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1));
                    rmlexp(iaddr(1,ibox)-1+tmpidx) = rmlexp(iaddr(1,ibox)-1+tmpidx) + reshape([real(cmlexpi(:)),imag(cmlexpi(:))]',[],1); 
                    %  rescale the multipole expansion, tmp is complex
                    cmlexpi = reshape(gboxmexp(:,i,jbox),nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
                    tmpi = zeros(size(cmlexpi));
                    tmp(:,:,:,ithd) = mpscale_mex(nd,nterms(ilev+1),cmlexpi,rscpow,tmpi);
                    %
                    %cc        process up down for current box
                    mexpupf = zeros(nd,nexptot); mexpdownf = zeros(nd,nexptot);
                    [mexpupf,mexpdownf] = mpoletoexp_mex(nd,tmp(:,:,:,ithd),nterms(ilev+2),nlams,nfourier,...
                            nexptot,mexpupf,mexpdownf,rlsc);
                    mexpf1(:,:,ithd) = mexpupf; mexpf2(:,:,ithd) = mexpdownf;
                    gboxwexp(:,:,1,i,ithd) = ftophys_mex(nd,mexpf1(:,:,ithd),nlams,rlams,nfourier,...
                            nphysical,nthmax,gboxwexp(:,:,1,i,ithd),fexpe,fexpo); % mexp is cumulative inside ftophys
                    gboxwexp(:,:,2,i,ithd) = ftophys_mex(nd,mexpf2(:,:,ithd),nlams,rlams,nfourier,...
                            nphysical,nthmax,gboxwexp(:,:,2,i,ithd),fexpe,fexpo);
                    [pgboxwexp(:,:,jbox,1),pgboxwexp(:,:,jbox,2)] = processgboxudexp_mex(nd,gboxwexp(:,:,1,i,ithd),gboxwexp(:,:,2,i,ithd),i,nexptotp,...
                            pgboxwexp(:,:,jbox,1),pgboxwexp(:,:,jbox,2),xshift,yshift,zshift);
                    %
                    %cc        process north-south for current box
                    cmptmpi = zeros(nd,(nterms(ilev+2)+1),(2*nterms(ilev+2)+1));
                    cmptmpi = rotztoy_mex(nd,nterms(ilev+2),tmp(:,:,:,ithd),cmptmpi,rdminus);
                    cmptmpi = cmptmpi(:);
                    mptmp(:,ithd) = reshape([real(cmptmpi) imag(cmptmpi)]',[],1); % 
                    mexpupf = zeros(nd,nexptot); mexpdownf = zeros(nd,nexptot);
                    [mexpupf,mexpdownf] = mpoletoexp_mex(nd,cmptmpi,nterms(ilev+2),nlams,nfourier,...
                            nexptot,mexpupf,mexpdownf,rlsc);
                    mexpf1(:,:,ithd) = mexpupf; mexpf2(:,:,ithd) = mexpdownf;
                    gboxwexp(:,:,3,i,ithd) = ftophys_mex(nd,mexpf1(:,:,ithd),nlams,rlams,nfourier,...
                            nphysical,nthmax,gboxwexp(:,:,3,i,ithd),fexpe,fexpo);
                    gboxwexp(:,:,4,i,ithd) = ftophys_mex(nd,mexpf2(:,:,ithd),nlams,rlams,nfourier,...
                            nphysical,nthmax,gboxwexp(:,:,4,i,ithd),fexpe,fexpo);
                    [pgboxwexp(:,:,jbox,3),pgboxwexp(:,:,jbox,4)] = processgboxudexp_mex(nd,gboxwexp(:,:,3,i,ithd),gboxwexp(:,:,4,i,ithd),i,nexptotp,...
                            pgboxwexp(:,:,jbox,3),pgboxwexp(:,:,jbox,4),xshift,yshift,zshift);
                    %
                    %cc        process east-west for current box  
                    cmptmpi = zeros(nd,(nterms(ilev+2)+1),(2*nterms(ilev+2)+1));
                    cmptmpi = rotztox_mex(nd,nterms(ilev+2),tmp(:,:,:,ithd),cmptmpi,rdplus);
                    cmptmpi = cmptmpi(:);
                    mptmp(:,ithd) = reshape([real(cmptmpi) imag(cmptmpi)]',[],1);
                    mexpupf = zeros(nd,nexptot); mexpdownf = zeros(nd,nexptot);
                    [mexpupf,mexpdownf] = mpoletoexp_mex(nd,cmptmpi,nterms(ilev+2),nlams,nfourier,...
                            nexptot,mexpupf,mexpdownf,rlsc);
                    mexpf1(:,:,ithd) = mexpupf; mexpf2(:,:,ithd) = mexpdownf;
                    gboxwexp(:,:,5,i,ithd) = ftophys_mex(nd,mexpf1(:,:,ithd),nlams,rlams,nfourier,...
                            nphysical,nthmax,gboxwexp(:,:,5,i,ithd),fexpe,fexpo);
                    gboxwexp(:,:,6,i,ithd) = ftophys_mex(nd,mexpf2(:,:,ithd),nlams,rlams,nfourier,...
                            nphysical,nthmax,gboxwexp(:,:,6,i,ithd),fexpe,fexpo);
                    [pgboxwexp(:,:,jbox,5),pgboxwexp(:,:,jbox,6)] = processgboxudexp_mex(nd,gboxwexp(:,:,5,i,ithd),gboxwexp(:,:,6,i,ithd),i,nexptotp,...
                            pgboxwexp(:,:,jbox,5),pgboxwexp(:,:,jbox,6),xshift,yshift,zshift);
                    
                  end
                end
              end
            end
         end
       end
%        rmpolesort_f=importdata('mps_data.dat'); % this checks out
%        mpolesort_f=importdata('mps_data.dat');
%        gboxmexp_f=importdata('mps_data.dat'); % this has 1e-17 difference?
%        rmlexp_f=importdata('mps_data.dat');
%        tmp_f=importdata('mps_data.dat');
%        pgboxwexp_f=importdata('mps_data.dat');

       clear gboxfl gboxsubcenters gboxwexp
       clear gboxind gboxisort_tmp gboxisort gboxmexp
%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 1 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc
       disp(['=== STEP 1 (shift mp) ===*']);     
%c     step 1, shift incoming multipole expansion to the center
%c     of each leaf-node box
       test = 0;
       for ilev=2:nlevels
         for ibox=laddr(1,ilev+1):laddr(2,ilev+1)
           istart = isrcse(1,ibox);
           iend = isrcse(2,ibox);
           npts = iend-istart+1;
           nchild = itree(ipointer(4)+ibox-1);
           if((npts>0)&&(nchild==0)&&(list4ct(ibox)==0))
             for i = istart:iend
               cmlexpi = zeros(nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
               cmlexpi = l3dmpmp_mex(nd,rmpolesort(i),cmpolesort(:,i),...
                   mpolesort(impolesort(i):(impolesort(i)+nd*(mtermssort(i)+1)*(2*mtermssort(i)+1)-1)),mtermssort(i),...
                   scales(ilev+1),treecenters(:,ibox),...
                   cmlexpi,nterms(ilev+1),dc,lca); % be careful with scales, nterms
               %%% here pay attention to the variable type, is it real or complex...
               tmpidx = 1:2*(nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1));
               rmlexp(iaddr(1,ibox)-1+tmpidx) = rmlexp(iaddr(1,ibox)-1+tmpidx) + reshape([real(cmlexpi(:)),imag(cmlexpi(:))]',[],1);
             end
           end
         end         
       end    
%        rmlexp_f=importdata('mps_data.dat');
%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 2 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc
       disp(['=== STEP 2 (merge mp) ===*']);  
       for ilev=nlevels-1:-1:0
         for ibox = laddr(1,ilev+1):laddr(2,ilev+1)
            for i=1:8
               jbox = itree(ipointer(5)+8*(ibox-1)+i-1);
               if(jbox>0) 
                  istart = isrcse(1,jbox);
                  iend = isrcse(2,jbox);
                  npts = iend-istart+1;
                  if(npts>0) 
                    %%% input cmlexpi needs special attention
                    rmlexpi = reshape(rmlexp(iaddr(1,jbox):(iaddr(1,jbox)+2*nd*(nterms(ilev+2)+1)*(2*nterms(ilev+2)+1)-1)),2,[]); % twice num of real ml
                    cmlexpi = reshape(rmlexpi(1,:)'+1i*rmlexpi(2,:)',nd,nterms(ilev+2)+1,2*nterms(ilev+2)+1);
                    cmlexpi2 = zeros(nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
                    cmlexpi2 = l3dmpmp_mex(nd,scales(ilev+2),treecenters(:,jbox),...
                        cmlexpi,nterms(ilev+2),...
                        scales(ilev+1),treecenters(:,ibox),...
                        cmlexpi2,nterms(ilev+1),dc,lca);
                    %%% here pay attention to the variable type, is it real or complex...
                    %%% add contribution instead of overwrite
                    tmpidx = 1:2*(nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1));
                    rmlexp(iaddr(1,ibox)-1+tmpidx) = rmlexp(iaddr(1,ibox)-1+tmpidx) + reshape([real(cmlexpi2(:)),imag(cmlexpi2(:))]',[],1);
                  end
               end
            end
         end   
       end
%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 3 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc
       disp(['=== STEP 3 (mp to loc+formta+mpeval) ===*']);  
       mexp = zeros(nd,nexptotp,nboxes,6);
 %cc   zero out mexp
       for k=1:6
         for i=1:nboxes
           for j=1:nexptotp
             for idim=1:nd
               mexp(idim,j,i,k) = 0.0d0;
             end
           end
         end
       end           
%c     init uall,dall,...,etc arrays
       uall=zeros(200,nthd); dall=zeros(200,nthd); nall=zeros(120,nthd); sall=zeros(120,nthd);
       eall=zeros(72,nthd); wall=zeros(72,nthd);
       u1234=zeros(36,nthd); d5678=zeros(36,nthd); 
       n1256=zeros(24,nthd); s3478=zeros(24,nthd);
       e1357=zeros(16,nthd); w2468=zeros(16,nthd); 
       n12=zeros(20,nthd); n56=zeros(20,nthd); s34=zeros(20,nthd); s78=zeros(20,nthd);
       e13=zeros(20,nthd); e57=zeros(20,nthd); w24=zeros(20,nthd); w68=zeros(20,nthd);
       e1=zeros(5,nthd); e3=zeros(5,nthd); e5=zeros(5,nthd); e7=zeros(5,nthd);
       w2=zeros(5,nthd); w4=zeros(5,nthd); w6=zeros(5,nthd); w8=zeros(5,nthd);
       iboxsubcenters=zeros(3,8,nthd);
       iboxfl=zeros(2,8,nthd);
%c     figure out allocations needed for iboxsrcind
%c     and so on
       nmaxt = 0;
       for ibox=1:nboxes % this is not used...
         if(nlist3(ibox)>0) 
           istart = isrcse(1,ibox);
           iend = isrcse(2,ibox);
           npts = iend-istart+1;
           if(npts>nmaxt) 
             nmaxt = npts;
           end
         end
       end
       iboxsrcind=zeros(nmaxt,nthd);
       iboxisort=zeros(nmaxt,nthd);
       iboxisort_tmp=zeros(nmaxt,nthd);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     for ilev = 2:nlevels
%      ilev = 2;
       iboxlexp=zeros(nd*(nterms(ilev)+1)*(2*nterms(ilev)+1),8,nthd);
       rscpow(1) = 1.0d0/boxsize(ilev+1);
       rtmp = scales(ilev+1)/boxsize(ilev+1);
       for i=1:nterms(ilev+1)
         rscpow(i+1) = rscpow(i)*rtmp;
       end
       testval = [];
       for ibox=laddr(1,ilev+1):laddr(2,ilev+1)
         ithd = 0;
         ithd = ithd + 1;
         istart = isrcse(1,ibox); 
         iend = isrcse(2,ibox);
         npts = iend-istart+1;
         if(npts>0)
%c         rescale the multipole expansion
           % tmp is complex, rmlexp is real, input cmlexpi needs special attention
           rmlexpi = reshape(rmlexp(iaddr(1,ibox):(iaddr(1,ibox)+2*nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1)-1)),2,[]); % twice num of real ml
           cmlexpi = reshape(rmlexpi(1,:)'+1i*rmlexpi(2,:)',nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
           tmpi = zeros(size(cmlexpi));
           tmp(:,:,:,ithd) = mpscale_mex(nd,nterms(ilev+1),cmlexpi,rscpow,tmpi);
%cc        process up down for current box
           mexpupf = zeros(nd,nexptot); mexpdownf = zeros(nd,nexptot);
           [mexpupf,mexpdownf] = mpoletoexp_mex(nd,tmp(:,:,:,ithd),nterms(ilev+1),nlams,nfourier,...
                   nexptot,mexpupf,mexpdownf,rlsc);
           mexpf1(:,:,ithd) = mexpupf; mexpf2(:,:,ithd) = mexpdownf;
           mexp(:,:,ibox,1) = ftophys_mex(nd,mexpf1(:,:,ithd),nlams,rlams,nfourier,...
                   nphysical,nthmax,mexp(:,:,ibox,1),fexpe,fexpo); % mexp is cumulative inside ftophys
           mexp(:,:,ibox,2) = ftophys_mex(nd,mexpf2(:,:,ithd),nlams,rlams,nfourier,...
                   nphysical,nthmax,mexp(:,:,ibox,2),fexpe,fexpo);
%cc        process north-south for current box, mptmp is real... rotztoy takes in complex
           cmptmpi = zeros(nd,(nterms(ilev+1)+1),(2*nterms(ilev+1)+1));
           cmptmpi = rotztoy_mex(nd,nterms(ilev+1),tmp(:,:,:,ithd),cmptmpi,rdminus);
           cmptmpi = cmptmpi(:);
           mptmp(:,ithd) = reshape([real(cmptmpi) imag(cmptmpi)]',[],1); % probably no use
           mexpupf = zeros(nd,nexptot); mexpdownf = zeros(nd,nexptot);
           [mexpupf,mexpdownf] = mpoletoexp_mex(nd,cmptmpi,nterms(ilev+1),nlams,nfourier,...
                   nexptot,mexpupf,mexpdownf,rlsc);
           mexpf1(:,:,ithd) = mexpupf; mexpf2(:,:,ithd) = mexpdownf;
           mexp(:,:,ibox,3) = ftophys_mex(nd,mexpf1(:,:,ithd),nlams,rlams,nfourier,...
                   nphysical,nthmax,mexp(:,:,ibox,3),fexpe,fexpo);
           mexp(:,:,ibox,4) = ftophys_mex(nd,mexpf2(:,:,ithd),nlams,rlams,nfourier,...
                   nphysical,nthmax,mexp(:,:,ibox,4),fexpe,fexpo);
%cc        process east-west for current box  
           cmptmpi = zeros(nd,(nterms(ilev+1)+1),(2*nterms(ilev+1)+1));
           cmptmpi = rotztox_mex(nd,nterms(ilev+1),tmp(:,:,:,ithd),cmptmpi,rdplus);
           cmptmpi = cmptmpi(:);
           mptmp(:,ithd) = reshape([real(cmptmpi) imag(cmptmpi)]',[],1);
           mexpupf = zeros(nd,nexptot); mexpdownf = zeros(nd,nexptot);
           [mexpupf,mexpdownf] = mpoletoexp_mex(nd,cmptmpi,nterms(ilev+1),nlams,nfourier,...
                   nexptot,mexpupf,mexpdownf,rlsc);
           mexpf1(:,:,ithd) = mexpupf; mexpf2(:,:,ithd) = mexpdownf;
           mexp(:,:,ibox,5) = ftophys_mex(nd,mexpf1(:,:,ithd),nlams,rlams,nfourier,...
                   nphysical,nthmax,mexp(:,:,ibox,5),fexpe,fexpo);
           mexp(:,:,ibox,6) = ftophys_mex(nd,mexpf2(:,:,ithd),nlams,rlams,nfourier,...
                   nphysical,nthmax,mexp(:,:,ibox,6),fexpe,fexpo);
         end
       end
%        rmlexp_f=importdata('mps_data.dat');
       mexpf12(:,:,1) = mexpf1(:,:,ithd);
       mexpf12(:,:,2) = mexpf2(:,:,ithd);

       nuall=0; ndall=0; nnall=0; nsall=0; neall=0; nwall=0; 
       nu1234=0; nd5678=0; nn1256=0; ns3478=0; ne1357=0; nw2468=0; 
       nn12=0; nn56=0; ns34=0; ns78=0; ne13=0; ne57=0; nw24=0; nw68=0; 
       ne1=0; ne3=0; ne5=0; ne7=0; nw2=0; nw4=0; nw6=0; nw8=0; 
       rscpow(1) = 1.0d0;
       rtmp = scales(ilev+1)/boxsize(ilev+1);
       for i=1:nterms(ilev+1)
         rscpow(i+1) = rscpow(i)*rtmp;
       end       
       for ibox = laddr(1,ilev):laddr(2,ilev)
         ithd = 0;
         ithd = ithd + 1;
         npts = 0;
         nchild = itree(ipointer(4)+ibox-1);
         istart = isrcse(1,ibox);
         iend = isrcse(2,ibox);
         npts = npts + iend-istart+1;
         if((npts>0)&&(nchild>0))
           % if(npts.gt.0.and.nchild.gt.0) then call getpwlistallprocessudnsewexp0
           nborsi = itree((ipointer(7)+mnbors*(ibox-1)):(ipointer(7)+mnbors*(ibox-1)-1+itree(ipointer(6)+ibox-1)));
           ichildi = itree(ipointer(5):(ipointer(5)-1+8*nboxes));
           
           rmlexp = getpwlistallprocessudnsewexp0_mex(ibox,boxsize(ilev+1),...
              nboxes,itree(ipointer(6)+ibox-1),...
              nborsi,...
              nchild,ichildi,treecenters,...
              isep,nd,ilev,scales(ilev+1),nterms(ilev+1),iaddr,...
              rmlexp,lmptot,...
              rlams,whts,nlams,nfourier,nphysical,...
              nthmax,nexptot,nexptotp,mexp,mexpf1(:,:,ithd),mexpf2(:,:,ithd),...
              mexpp1(:,:,ithd),mexpp2(:,:,ithd),...
              mexppall(:,:,:,ithd),rdplus,rdminus,...
              xshift,yshift,zshift,fexpback,nn,rlsc,rscpow,...
              pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4); 
         end
         if((nlist3(ibox)>0)&&(npts>0)) 
           disp(['this is actually processed...']);  
           for i=1:8
              iboxlexp(:,i,ithd) = 0;
           end
           iboxlexp(:,:,ithd) = getlist3pwallprocessudnsewexp0_mex(ibox,boxsize(ilev+1),...
                 nboxes,nlist3(ibox),list3(1:nlist3(ibox),ibox),...
                 isep,treecenters,nd,...
                 nterms(ilev+1),iboxlexp(:,:,ithd),rlams,...
                 whts,nlams,nfourier,nphysical,nthmax,...
                 nexptot,nexptotp,mexp,...
                 mexpf1(:,:,ithd),mexpf2(:,:,ithd),...
                 mexpp1(:,:,ithd),mexpp2(:,:,ithd),...
                 mexppall(:,:,1,ithd),mexppall(:,:,2,ithd),...
                 xshift,yshift,zshift,fexpback,nn,rlsc,rscpow,...
                 rdplus,rdminus);
           
           istart = isrcse(1,ibox);
           iend = isrcse(2,ibox);
           npts = iend-istart+1;
           if(npts>0) 
             iboxsrcindi = iboxsrcind(:,ithd);
             iboxfli = iboxfl(:,:,ithd);
             iboxsubcentersi = iboxsubcenters(:,:,ithd);
             [iboxsrcind(:,ithd),iboxfl(:,:,ithd),iboxsubcenters(:,:,ithd)] = ...
                  subdividebox_mex(cmpolesort(:,istart),npts,treecenters(:,ibox),boxsize(ilev+1),iboxsrcindi,iboxfli,iboxsubcentersi);
             for i=istart:iend
                 iboxisort_tmp(i-istart+1,ithd) = i;
             end
             iboxisort_tmpi = iboxisort_tmp(:,ithd);
             iboxisorti = iboxisort(:,ithd);
             iboxsrcindi = iboxsrcind(:,ithd);
             iboxisort(:,ithd) = ireorderf_mex(1,npts,iboxisort_tmpi',iboxisorti',iboxsrcindi)';
             for i=1:8
               if(iboxfl(1,i,ithd)>0) 
                 jstart=iboxfl(1,i,ithd);
                 jend=iboxfl(2,i,ithd);
                 npts0=jend-jstart+1;
                 for j = jstart:jend
                   k = iboxisort(j,ithd);
                   localsortj = localsort(impolesort(k):(impolesort(k)-1+nd*(mtermssort(k)+1)*(2*mtermssort(k)+1)));
                   localsortj = l3dlocloc_mex(nd,scales(ilev+1),...
                         iboxsubcenters(:,i,ithd),iboxlexp(:,i,ithd),...
                         nterms(ilev+1),...
                         rmpolesort(k),cmpolesort(:,k),...
                         localsortj, mtermssort(k),...
                         dc,lca);
                   localsort(impolesort(k):(impolesort(k)-1+nd*(mtermssort(k)+1)*(2*mtermssort(k)+1))) = localsortj(:);
                 end
               end
             end
                
%              keyboard
           end
         end
       end

%        
%        iboxlexp_f=importdata('mps_data.dat');

     end
%     localsort_f=importdata('mps_data.dat'); % localsort check...
%     rmlexp_f=importdata('mps_data.dat'); % rmlexp check...
%
%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 4 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc     
     disp(['=== Step 4 (split loc) ===*']);    
     for ilev = 2:nlevels-1
       for ibox = laddr(1,ilev+1):laddr(2,ilev+1)
         istart = isrcse(1,ibox);
         iend = isrcse(2,ibox);
         npts = iend-istart+1;
         if(npts>0) 
           for i=1:8
             jbox = itree(ipointer(5)+8*(ibox-1)+i-1);
             if(jbox>0) % pay attention to real or complex
               rmlexpi = reshape(rmlexp(iaddr(2,ibox):(iaddr(2,ibox)-1+2*nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1))),2,[]);
               cmlexpi = reshape(rmlexpi(1,:)'+1i*rmlexpi(2,:)',nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
               % assign 0 value to rmlexpj, then add output from l3dlocloc back to rmlexp
               rmlexpj = 0*reshape(rmlexp(iaddr(2,jbox):(iaddr(2,jbox)-1+2*nd*(nterms(ilev+2)+1)*(2*nterms(ilev+2)+1))),2,[]);
               cmlexpj = reshape(rmlexpj(1,:)'+1i*rmlexpj(2,:)',nd,nterms(ilev+2)+1,2*nterms(ilev+2)+1);
               cmlexpj = l3dlocloc_mex(nd,scales(ilev+1),...
                     treecenters(:,ibox),cmlexpi,...
                     nterms(ilev+1),scales(ilev+2),treecenters(:,jbox),...
                     cmlexpj,nterms(ilev+2),dc,lca);
               rmlexpj = [real(cmlexpj(:)) imag(cmlexpj(:))]';
               rmlexp(iaddr(2,jbox):(iaddr(2,jbox)-1+2*nd*(nterms(ilev+2)+1)*(2*nterms(ilev+2)+1))) = rmlexp(iaddr(2,jbox):(iaddr(2,jbox)-1+2*nd*(nterms(ilev+2)+1)*(2*nterms(ilev+2)+1))) + rmlexpj(:);
             end
           end
         end
       end
     end
%      rmlexp_f=importdata('mps_data.dat'); % rmlexp check...
%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 5 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc
% c      ... step 5, shift leaf box loc exp to mpole center
     disp(['=== Step 5 (LOC to CEN) ===*']);
     for ilev = 0:nlevels
       for ibox = laddr(1,ilev+1):laddr(2,ilev+1)
         nchild = itree(ipointer(4)+ibox-1);
         if(nchild==0) % leaf box loc exp to mpole center
           istart = isrcse(1,ibox);
           iend = isrcse(2,ibox);
           npts = iend - istart + 1;
           for i = istart:iend % pay attention to real or complex...
             rmlexpi = reshape(rmlexp(iaddr(2,ibox):(iaddr(2,ibox)-1+2*nd*(nterms(ilev+1)+1)*(2*nterms(ilev+1)+1))),2,[]);
             cmlexpi = reshape(rmlexpi(1,:)'+1i*rmlexpi(2,:)',nd,nterms(ilev+1)+1,2*nterms(ilev+1)+1);
             % assign 0 value to localsorti, localsort is complex already
             localsorti = 0*localsort(impolesort(i):(impolesort(i)-1+nd*(mtermssort(i)+1)*(2*mtermssort(i)+1)));
             localsorti = l3dlocloc_mex(nd, scales(ilev+1),...
                   treecenters(:,ibox), cmlexpi,...
                   nterms(ilev+1), rmpolesort(i), cmpolesort(:,i),...
                   localsorti, mtermssort(i),...
                   dc,lca);
             localsort(impolesort(i):(impolesort(i)-1+nd*(mtermssort(i)+1)*(2*mtermssort(i)+1))) = localsort(impolesort(i):(impolesort(i)-1+nd*(mtermssort(i)+1)*(2*mtermssort(i)+1))) + localsorti(:);
           end
         end
       end
     end
%      localsort_f = importdata('mps_data.dat');

%cccccc
%cccccc used to be insdie lfmm3dmain_mps       
%cccccc STEP 6 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc   
%cc     directly evaluate potential at sources and targets
     disp(['=== STEP 6 (direct) =====*']);
     for ilev=0:nlevels
       for ibox = laddr(1,ilev+1):laddr(2,ilev+1)
         istart = isrcse(1,ibox);
         iend = isrcse(2,ibox);
         npts0 = iend-istart+1;
         for iloc = istart:iend
           for i = 1:nlist1(ibox)
             jbox = list1(i,ibox);
             jstart = isrcse(1,jbox);
             jend = isrcse(2,jbox);
             npts = jend - jstart+1;
             for j = jstart:jend
               d = (cmpolesort(1,j)-cmpolesort(1,iloc))^2 + (cmpolesort(2,j)-cmpolesort(2,iloc))^2 + (cmpolesort(3,j)-cmpolesort(3,iloc))^2;
               d = sqrt(d);
               if (d > thresh)
                 mpolesortj = mpolesort(impolesort(j):(impolesort(j)-1+nd*(mtermssort(j)+1)*(2*mtermssort(j)+1)));
                 % assign 0 value to localsortj
                 localsortj = 0*localsort(impolesort(iloc):(impolesort(iloc)-1+nd*(mtermssort(iloc)+1)*(2*mtermssort(iloc)+1)));
                 localsortj = l3dmploc_mex(nd, rmpolesort(j),...
                       cmpolesort(:,j),...
                       mpolesortj, mtermssort(j),...
                       rmpolesort(iloc), cmpolesort(:,iloc),...
                       localsortj,...
                       mtermssort(iloc),dc,lca);
                 localsort(impolesort(iloc):(impolesort(iloc)-1+nd*(mtermssort(iloc)+1)*(2*mtermssort(iloc)+1))) = localsort(impolesort(iloc):(impolesort(iloc)-1+nd*(mtermssort(iloc)+1)*(2*mtermssort(iloc)+1))) + localsortj(:);
               end
             end
           end
         end
       end
     end
% localsort_f = importdata('mps_data.dat');

%cccccc       
%cccccc used to be insdie lfmm3d_mps, after lfmm3dmain_mps
%cccccc sor local expansion ccccccccccccccccccccccccccccccccccccccccccccccc
%cccccc
     for i = 1:nmpole
       mt = mtermssort(i);
       ilen = (mt+1)*(2*mt+1);
       ijk = 1;
       for j = 1:ilen
         for l = 1:nd
           local(impole(isrc(i))+ijk-1)=localsort(impolesort(i)+ijk-1);
           ijk = ijk+1;
         end
       end
     end
% local_f = importdata('mps_data.dat');

%%%%%% post process       
       interms = (mterms(1)+1)*(2*mterms(1)+1);
       ilocal = zeros(mterms(1)+1,2*mterms(1)+1);
       ilocal = 0;
       tmp_vec = zeros(interms,1);
       pot2 = zeros(nd,nmpole);
       npts = 1;
       for i = 1:nmpole
         pot2tmp = l3dtaevalp_mex(nd, rmpole(i),cmpole(:,i), local(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1), mterms(i), source(:,i), npts, pot2(:,i), wlege, nlege);
         pot2(:,i) = pot2tmp;
       end

       abs(pot(1,1:10)-pot2(1,1:10))
       
       keyboard

