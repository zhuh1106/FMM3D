function local = l3dmploc_mex(nd,sc1,x0y0z0,mpole,nterms,sc2,xnynzn,local,nterms2,dc,lca)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
ndnterms2 = nd*(nterms2+1);
nterms2t2p1 = 2*nterms2+1;
local = reshape(local,ndnterms2,nterms2t2p1);
lca1 = lca+1;
mex_id_ = 'l3dmploc(i int64_t[x], i double[x], i double[x], i dcomplex[xx], i int64_t[x], i double[x], i double[x], io dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x])';
[local] = laprouts3d(mex_id_, nd, sc1, x0y0z0, mpole, nterms, sc2, xnynzn, local, nterms2, dc, lca, 1, 1, 3, ndnterms, ntermst2p1, 1, 1, 3, ndnterms2, nterms2t2p1, 1, lca1, lca1, 1);
local = reshape(local,nd,nterms2+1,nterms2t2p1);
end

