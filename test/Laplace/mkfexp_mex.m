function [fexpe,fexpo,fexpback] = mkfexp_mex(nlambs,numfour,numphys,fexpe,fexpo,fexpback)
nn = numel(fexpe);
mex_id_ = 'mkfexp(i int[x], i int[x], i int[x], io dcomplex[x], io dcomplex[x], io dcomplex[x])';
[fexpe, fexpo, fexpback] = laprouts3d(mex_id_, nlambs, numfour, numphys, fexpe, fexpo, fexpback, 1, nlambs, nlambs, nn, nn, nn);
end