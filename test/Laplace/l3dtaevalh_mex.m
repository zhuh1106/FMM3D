function [pot,grad,hess] = l3dtaevalh_mex(nd,rscale,center,local,nterms,ztarg,ntarg,pot,grad,hess,scarray)
nd3 = 3*nd;
nd6 = 6*nd;
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
local = reshape(local,ndnterms,ntermst2p1);
grad = reshape(grad,nd3,ntarg);
hess = reshape(hess,nd6,ntarg);
nscarray = length(scarray);
mex_id_ = 'l3dtaevalh(i int[x], i double[x], i double[x], i dcomplex[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pot, grad, hess] = laprouts3d(mex_id_, nd, rscale, center, local, nterms, ztarg, ntarg, pot, grad, hess, scarray, 1, 1, 3, ndnterms, ntermst2p1, 1, 3, ntarg, 1, nd, ntarg, nd3, ntarg, nd6, ntarg, nscarray);
grad = reshape(grad,nd,3,ntarg);
hess = reshape(hess,nd,6,ntarg);
end

