%
%
%
      addpath('../../matlab/')
      clear

      done=1;
      pi=4.0*atan(done);
      
%       set scale parameter neq 1 to make sure it is used correctly.

      bsize = 1.0d0;
      rscale1 = bsize/4;
      rscale2 = bsize/2;

      rscale1 = 0.55;
      rscale2 = rscale1*2;

      nd = 1;

%       set center for original multipole expansion

      c0(1)=0.248d0;
      c0(2)=0.251d0;
      c0(3)=0.249d0;

%       create two charge sources.

      sources(1,1)=c0(1)+0.24d0;
      sources(2,1)=c0(2)+ 0.25d0;
      sources(3,1)=c0(3)+ 0.25d0;
      for idim=1:nd
        charge(idim,1) = 1.0d0;
        dipvec(idim,1,1) = 1.0d0/2;
        dipvec(idim,2,1) = 2.0d0/3;
        dipvec(idim,3,1) = 3.0d0/4;
      end


      ns = 1;
      sources(1,2)=c0(1)+0.25d0;
      sources(2,2)=c0(2)-0.25d0;
      sources(3,2)=c0(3)-0.25d0;
      for idim=1:nd
        charge(idim,2)= -1.0d0;  
        dipvec(idim,1,2) = 3.0d0/5;
        dipvec(idim,2,2) = 2.0d0/6;
        dipvec(idim,3,2) = 1.0d0/7;
      end
      ns = 2;

%       create center for shifted expansions
%       (location jiggled for good measure)
%
      c1(1)= 0.015d0;
      c1(2)= 0.012d0;
      c1(3)= 0.013d0;
      c2(1)= 2.015d0;
      c2(2)= 0.012d0;
      c2(3)= 0.013d0;
      c3(1)= c2(1)-0.25d0;
      c3(2)= c2(2)-0.249d0;
      c3(3)= c2(3)-0.251d0;
%
%       create target
%
      ztrg(1)=c3(1)-0.245d0;
      ztrg(2)=c3(2)-0.25d0;
      ztrg(3)=c3(3)-0.25d0;
      ztrgs(1,1)= ztrg(1);
      ztrgs(2,1)= ztrg(2);
      ztrgs(3,1)= ztrg(3);
      ztrgs(1,2)= ztrg(1)+0.05d0;
      ztrgs(2,2)= ztrg(2)+0.05d0;
      ztrgs(3,2)= ztrg(3)+0.05d0;
      nt = 2;

     
      ntest = 5;
      for i=1:ntest
        ipass(i) = 0; 
      end

%
%       direct calculation:
%
      for i=1:nt
        for idim=1:nd
          opots(idim,i) = 0;
          oflds(idim,1,i) = 0;
          oflds(idim,2,i) = 0;
          oflds(idim,3,i) = 0;
          ohesss(idim,1,i) = 0;
          ohesss(idim,2,i) = 0;
          ohesss(idim,3,i) = 0;
          ohesss(idim,4,i) = 0;
          ohesss(idim,5,i) = 0;
          ohesss(idim,6,i) = 0;
        end
      end

     thresh = 1.0d-15;

     [opots, oflds, ohesss] = l3ddirectcdh_mex(nd, sources, charge, dipvec, ns, ztrgs, nt, opots, oflds, ohesss, thresh);


    % copied from lfmm3d test
    srcinfo.sources = sources;
    srcinfo.charges = charge;
    srcinfo.dipoles = dipvec;
    srcinfo.nd = nd;
    ttmp = ztrgs;
    U2 = l3ddir(srcinfo,ttmp,3);
    U2.pot = U2.pottarg;
    max(abs(U2.pot(:) - opots(:)))

      eps = 0.5d-12;
      nterms = 0; nterms2 = 0; nterms3 = 0;
      nterms = l3dterms_mex(eps, nterms);
      nterms2 = l3dterms_mex(eps, nterms2);
      nterms3 = l3dterms_mex(eps, nterms3);

      rconv1 = 1.0d0/sqrt(3.0d0);
      rconv2 = sqrt(3.0d0)/2.0d0;
      rconv3 = 0.75d0;

      nlege = 100;
      lw7 = 100000;
      wlege = zeros(2*(nlege+1)^2,1);
      lused7 = 0;
      [wlege,lused7]=ylgndrfwini_mex(nlege,wlege,lw7,lused7);
      scarray_mp=zeros(9*(nterms+1)^2,1); % why 10 in fortran subroutine?
      scarray_mp=l3dmpevalhessdini_mex(nterms,scarray_mp);
      iuse = 9*(nterms+1)^2-3*1-5*(1+3); % double check... 
      scarray_loc = zeros(iuse,1);
      scarray_loc = l3dtaevalhessdini_mex(nterms,scarray_loc);
      % l3dtaevalhessdini(nterms,scarray_loc)

%
%
%       create multipole expansion:
%
      mpole1 = zeros(nd,nterms+1,2*nterms+1);
      mpole1 = mpzero_mex(nd,mpole1,nterms);
      mpole1 = l3dformmpcd_mex(nd,rscale1,sources,charge,dipvec,ns,c0,nterms,mpole1,wlege,nlege);

      
      for i=1:nt
        for idim=1:nd
          pots(idim,i) = 0;
          flds(idim,1,i) = 0;
          flds(idim,2,i) = 0;
          flds(idim,3,i) = 0;
          hesss(idim,1,i) = 0;
          hesss(idim,2,i) = 0;
          hesss(idim,3,i) = 0;
          hesss(idim,4,i) = 0;
          hesss(idim,5,i) = 0;
          hesss(idim,6,i) = 0;
        end
      end 

     [pots,flds,hesss] = l3dmpevalh_mex(nd,rscale1,c0,mpole1,nterms,ztrgs,nt,pots,flds,hesss,thresh,scarray_mp);    

      pots - opots
%
%    mpmp shift
%
      nn = nterms;
      nn = max(nn,nterms2);
      nn = max(nn,nterms3);
      nn = 2*nn + 10;
      dc = zeros(nn+1,nn+1);
      dc = getsqrtbinomialcoeffs_mex(nn,dc);

      mpole2 = zeros(nd,nterms2+1,2*nterms2+1);
      mpole2 = mpzero_mex(nd,mpole2,nterms2);
      mpole2 = l3dmpmp_mex(nd,rscale1,c0,mpole1,nterms,rscale2,c1,mpole2,nterms2,dc,nn);

      for i=1:nt
        for idim=1:nd
          pots(idim,i) = 0;
          flds(idim,1,i) = 0;
          flds(idim,2,i) = 0;
          flds(idim,3,i) = 0;
          hesss(idim,1,i) = 0;
          hesss(idim,2,i) = 0;
          hesss(idim,3,i) = 0;
          hesss(idim,4,i) = 0;
          hesss(idim,5,i) = 0;
          hesss(idim,6,i) = 0;
        end
      end

      [pots,flds,hesss] = l3dmpevalh_mex(nd,rscale2,c1,mpole2,nterms,ztrgs,nt,pots,flds,hesss,thresh,scarray_mp);

      pots - opots
      flds - oflds
      hesss - ohesss

%
%     convert to local
%
      locexp2 = zeros(nd,nterms2+1,2*nterms2+1);
      locexp2 = mpzero_mex(nd,locexp2,nterms2);
      locexp2 = l3dmploc_mex(nd,rscale2,c1,mpole2,nterms2,rscale2,c2,locexp2,nterms2,dc,nn);

      for i=1:nt
        for idim=1:nd
          pots(idim,i) = 0;
          flds(idim,1,i) = 0;
          flds(idim,2,i) = 0;
          flds(idim,3,i) = 0;
          hesss(idim,1,i) = 0;
          hesss(idim,2,i) = 0;
          hesss(idim,3,i) = 0;
          hesss(idim,4,i) = 0;
          hesss(idim,5,i) = 0;
          hesss(idim,6,i) = 0;
        end
      end

      [pots,flds,hesss] = l3dtaevalh_mex(nd,rscale2,c2,locexp2,nterms2,ztrgs,nt,pots,flds,hesss,scarray_loc);
      pots-opots
      flds-oflds
      hesss-ohesss

%
%     shift local 

      locexp1 =zeros(nd,nterms3+1,2*nterms3+1);
      locexp1 = mpzero_mex(nd,locexp1,nterms3);

      locexp1 = l3dlocloc_mex(nd,rscale2,c2,locexp2,nterms2,rscale1,c3,locexp1,nterms3,dc,nn);

      for i=1:nt
        for idim=1:nd
          pots(idim,i) = 0;
          flds(idim,1,i) = 0;
          flds(idim,2,i) = 0;
          flds(idim,3,i) = 0;
          hesss(idim,1,i) = 0;
          hesss(idim,2,i) = 0;
          hesss(idim,3,i) = 0;
          hesss(idim,4,i) = 0;
          hesss(idim,5,i) = 0;
          hesss(idim,6,i) = 0;
        end
      end  

      [pots,flds,hesss] = l3dtaevalh_mex(nd,rscale1,c3,locexp1,nterms3,ztrgs,nt,pots,flds,hesss,scarray_loc);
      pots-opots
      flds-oflds
      hesss-ohesss

%
%
%
%    create local exp from sources 
      locexp1 = mpzero_mex(nd,locexp1,nterms3);
      locexp1 = l3dformtacd_mex(nd,rscale1,sources,charge,dipvec,ns,c3,nterms3,locexp1,wlege,nlege);

      for i=1:nt
        for idim=1:nd
          pots(idim,i) = 0;
          flds(idim,1,i) = 0;
          flds(idim,2,i) = 0;
          flds(idim,3,i) = 0;
          hesss(idim,1,i) = 0;
          hesss(idim,2,i) = 0;
          hesss(idim,3,i) = 0;
          hesss(idim,4,i) = 0;
          hesss(idim,5,i) = 0;
          hesss(idim,6,i) = 0;
        end
      end

      [pots,flds,hesss] = l3dtaevalh_mex(nd,rscale1,c3,locexp1,nterms3,ztrgs,nt,pots,flds,hesss,scarray_loc);

      pots-opots
      flds-oflds
      hesss-ohesss
      
      keyboard
