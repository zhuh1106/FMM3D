function mpole = mpzero_mex(nd,mpole,nterms)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
mex_id_ = 'mpzero(i int64_t[x], io dcomplex[xx], i int64_t[x])';
[mpole] = laprouts3d(mex_id_, nd, mpole, nterms, 1, ndnterms, ntermst2p1, 1);
mpole = reshape(mpole,nd,nterms+1,ntermst2p1);
end

