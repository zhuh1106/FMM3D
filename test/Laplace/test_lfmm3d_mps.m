% test lfmm3d_mps
%
%
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
      
      
       n1 = 4;
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
           charge(idim,i) = hkrand_mex(0);
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
      
       ntm = 7;
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
         mpoletmp =  l3dformmpc_mex(nd, rscale, source(:,i), charge(:,i), ns1, cmpole(:,i), mterms(i), mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1),wlege, nlege);
         mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) = mpoletmp(:);
%          keyboard
       end

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
       local = lfmm3d_mps0_mex(nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);
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


       
       keyboard


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

