function [xs,ys,zs] = mkexps_mex(rlams,nlambs,numphys,nexptotp,xs,ys,zs)
mex_id_ = 'mkexps(i double[x], i int[x], i int[x], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx])';
[xs, ys, zs] = laprouts3d(mex_id_, rlams, nlambs, numphys, nexptotp, xs, ys, zs, nlambs, 1, nlambs, 1, 11, nexptotp, 11, nexptotp, 5, nexptotp);
end

