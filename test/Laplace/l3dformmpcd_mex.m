function mpole = l3dformmpcd_mex(nd,rscale,sources,charge,dipvec,ns,center,nterms,mpole,wlege,nlege)
nd3 = 3*nd;
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
% nlege1 = nlege+1; % although wlege is defined with dimension nlege1^2, it is actually 2*(nlege+1)^2 
nwlege = numel(wlege);
dipvec = reshape(dipvec,nd3,ns);
mpole = reshape(mpole,ndnterms,ntermst2p1);
mex_id_ = 'l3dformmpcd(i int64_t[x], i double[x], i double[xx], i double[xx], i double[xx], i int64_t[x], i double[x], i int64_t[x], io dcomplex[xx], i double[x], i int64_t[x])';
[mpole] = laprouts3d(mex_id_, nd, rscale, sources, charge, dipvec, ns, center, nterms, mpole, wlege, nlege, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, 1, ndnterms, ntermst2p1, nwlege, 1);
mpole = reshape(mpole,nd,nterms+1,2*nterms+1);
end

