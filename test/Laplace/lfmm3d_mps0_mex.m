function local = lfmm3d_mps0_mex(nd, eps ,nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier)
ntmp = numel(mpole);
mtmp = numel(local);
mex_id_ = 'lfmm3dmps0(i int[x], i double[x], i int[x], i double[xx], i double[x], i int[x], i dcomplex[x], i int[x], io dcomplex[x], i int[x])';
[local] = laprouts3d(mex_id_, nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier, 1, 1, 1, 3, nmpole, nmpole, nmpole, ntmp, nmpole, mtmp, 1);
end

