function pot = l3dtaevalp_mex(nd,rscale,center,mpole,nterms,ztarg,ntarg,pot,wlege,nlege)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
nwlege = numel(wlege);
mex_id_ = 'l3dtaevalp(i int64_t[x], i double[x], i double[x], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io double[xx], i double[x], i int64_t[x])';
[pot] = laprouts3d(mex_id_, nd, rscale, center, mpole, nterms, ztarg, ntarg, pot, wlege, nlege, 1, 1, 3, ndnterms, ntermst2p1, 1, 3, ntarg, 1, nd, ntarg, nwlege, 1);
end

