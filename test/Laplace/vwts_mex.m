function [x,w] = vwts_mex(x,w,n)
mex_id_ = 'vwts(io double[x], io double[x], i int[x])';
[x, w] = laprouts3d(mex_id_, x, w, n, n, n, 1);
end

