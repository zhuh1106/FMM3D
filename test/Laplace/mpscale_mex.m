function mpoleout = mpscale_mex(nd,nterms,mpolein,rsc,mpoleout)
nterms1 = nterms+1;
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpolein = reshape(mpolein,ndnterms,ntermst2p1);
mpoleout = reshape(mpoleout,ndnterms,ntermst2p1);
mex_id_ = 'mpscale(i int64_t[x], i int64_t[x], i dcomplex[xx], i double[x], io dcomplex[xx])';
[mpoleout] = laprouts3d(mex_id_, nd, nterms, mpolein, rsc, mpoleout, 1, 1, ndnterms, ntermst2p1, nterms1, ndnterms, ntermst2p1);
mpoleout = reshape(mpoleout,nd,nterms+1,2*nterms+1);
end

