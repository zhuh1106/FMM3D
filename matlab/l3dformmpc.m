function mpole = l3dformmpc(nd,rscale,sources,charge,ns,center,nterms,mpole,wlege,nlege)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
nwlege = numel(wlege);
mex_id_ = 'l3dformmpc(i int[x], i double[x], i double[xx], i double[xx], i int[x], i double[x], i int[x], io dcomplex[xx], i double[x], i int[x])';
[mpole] = fmm3d(mex_id_, nd, rscale, sources, charge, ns, center, nterms, mpole, wlege, nlege, 1, 1, 3, ns, nd, ns, 1, 3, 1, ndnterms, ntermst2p1, nwlege, 1);
mpole = reshape(mpole,nd,nterms+1,2*nterms+1);
end

% ---------------------------------------------------------------------
