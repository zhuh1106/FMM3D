function mexpphys = ftophys_mex(nd,mexpf,nlambs,rlams,numfour,numphys,nthmax,mexpphys,fexpe,fexpo)
nexptot = size(mexpf,2);
nn = numel(fexpe);
nexptotp = size(mexpphys,2);
mex_id_ = 'ftophys(i int64_t[x], i dcomplex[xx], i int64_t[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], io dcomplex[xx], i dcomplex[x], i dcomplex[x])';
[mexpphys] = laprouts3d(mex_id_, nd, mexpf, nlambs, rlams, numfour, numphys, nthmax, mexpphys, fexpe, fexpo, 1, nd, nexptot, 1, nlambs, nlambs, nlambs, 1, nd, nexptotp, nn, nn);
end

