function [w,lused] = ylgndrfwini_mex(nmax, w, lw, lused)
nw = length(w);
mex_id_ = 'ylgndrfwini(i int[x], io double[x], i int[x], io int[x])';
[w, lused] = laprouts3d(mex_id_, nmax, w, lw, lused, 1, nw, 1, 1);
end

