function mpole = h3dformmpc(nd, zk, rscale, sources, charge, ns, center, nterms, wlege, nlege)
%
%
mpole_dim = nd*(nterms+1)*(2*nterms+1);
s_dim = numel(sources);
c_dim = nd*s_dim/3;
mpole = zeros(mpole_dim,1);
mex_id_ = 'h3dformmpc(i int[x], i dcomplex[x], i double[x], i double[x], i dcomplex[x], i int[x], i double[x], i int[x], io dcomplex[x], i double[], i int[x])';
[mpole] = fmps(mex_id_, nd, zk, rscale, sources, charge, ns, center, nterms, mpole, wlege, nlege, 1, 1, 1, s_dim, c_dim, 1, 3, 1, mpole_dim, 1);

