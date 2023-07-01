function mpole = l3dformmpc_mex(nd,rscale,sources,charge,ns,center,nterms,mpole,wlege,nlege)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
nwlege = numel(wlege);
mex_id_ = 'l3dformmpc(i int[x], i double[x], i double[xx], i double[xx], i int[x], i double[x], i int[x], io dcomplex[xx], i double[x], i int[x])';
[mpole] = laprouts3d(mex_id_, nd, rscale, sources, charge, ns, center, nterms, mpole, wlege, nlege, 1, 1, 3, ns, nd, ns, 1, 3, 1, ndnterms, ntermst2p1, nwlege, 1);
mpole = reshape(mpole,nd,nterms+1,2*nterms+1);
end

%@function [pot,grad,hess,pottarg,gradtarg,hesstarg] = lfmm3d_mex(nd,eps,ns,source,ifcharge,charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,nt,targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
%% need to take care of ns or nt = 0 ?
%nd3 = nd*3;
%nd6 = nd*6;
%dipvec = reshape(dipvec,nd3,ns);
%grad = reshape(grad,nd3,ns);
%hess = reshape(hess,nd6,ns);
%gradtarg = reshape(gradtarg,nd3,nt);
%hesstarg = reshape(hesstarg,nd6,nt);
%# FORTRAN lfmm3d(int[1] nd, double[1] eps, int[1] ns, double[3,ns] source, int[1] ifcharge, double[nd,ns] charge, int[1] ifdipole, double[nd3,ns] dipvec, int[1] iper, int[1] ifpgh, inout double[nd,ns] pot, inout double[nd3,ns] grad, inout double[nd6,ns] hess, int[1] nt, double[3,nt] targ, int[1] ifpghtarg, inout double[nd,nt] pottarg, inout double[nd3,nt] gradtarg, inout double[nd6,nt] hesstarg, int[1] ier);
%grad = reshape(grad,nd,3,ns);
%hess = reshape(hess,nd,6,ns);
%gradtarg = reshape(gradtarg,nd,3,nt);
%hesstarg = reshape(hesstarg,nd,6,nt);
%end

