function [pot,grad,hess] = l3dmpevalh_mex(nd,rscale,center,mpole,nterms,ztarg,nt,pot,grad,hess,thresh,scarray)
nd3 = 3*nd;
nd6 = 6*nd;
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
grad = reshape(grad,nd3,nt);
hess = reshape(hess,nd6,nt);
nscarray = length(scarray);
mex_id_ = 'l3dmpevalh(i int64_t[x], i double[x], i double[x], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io double[xx], io double[xx], io double[xx], i double[x], i double[x])';
[pot, grad, hess] = laprouts3d(mex_id_, nd, rscale, center, mpole, nterms, ztarg, nt, pot, grad, hess, thresh, scarray, 1, 1, 3, ndnterms, ntermst2p1, 1, 3, nt, 1, nd, nt, nd3, nt, nd6, nt, 1, nscarray);
grad = reshape(grad,nd,3,nt);
hess = reshape(hess,nd,6,nt);
end

