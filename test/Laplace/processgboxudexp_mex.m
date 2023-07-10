function [mexpuall,mexpdall] = processgboxudexp_mex(nd,mexpugbox,mexpdgbox,jbox,nexptotp,mexpuall,mexpdall,xs,ys,zs)
mex_id_ = 'processgboxudexp(i int64_t[x], i dcomplex[xx], i dcomplex[xx], i int64_t[x], i int64_t[x], io dcomplex[xx], io dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i double[xx])';
[mexpuall, mexpdall] = laprouts3d(mex_id_, nd, mexpugbox, mexpdgbox, jbox, nexptotp, mexpuall, mexpdall, xs, ys, zs, 1, nd, nexptotp, nd, nexptotp, 1, 1, nd, nexptotp, nd, nexptotp, 11, nexptotp, 11, nexptotp, 5, nexptotp);
end

