function local = l3dlocloc_mex(nd,sc1,x0y0z0,locold,nterms,sc2,xnynzn,local,nterms2,dc,lda)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
locold = reshape(locold,ndnterms,ntermst2p1);
ndnterms2 = nd*(nterms2+1);
nterms2t2p1 = 2*nterms2+1;
local = reshape(local,ndnterms2,nterms2t2p1);
lda1 = lda+1;
mex_id_ = 'l3dlocloc(i int[x], i double[x], i double[x], i dcomplex[xx], i int[x], i double[x], i double[x], io dcomplex[xx], i int[x], i double[xx], i int[x])';
[local] = laprouts3d(mex_id_, nd, sc1, x0y0z0, locold, nterms, sc2, xnynzn, local, nterms2, dc, lda, 1, 1, 3, ndnterms, ntermst2p1, 1, 1, 3, ndnterms2, nterms2t2p1, 1, lda1, lda1, 1);
local = reshape(local,nd,nterms2+1,nterms2t2p1);
end

