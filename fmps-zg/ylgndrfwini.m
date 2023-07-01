function [wlege, lused] = ylgndrfwini(nlege,wlege,lw,lused)
%
%
% have to use inout... maybe modify this fortran routine, or is this standard practice... Hai
mex_id_ = 'ylgndrfwini(i int[x], io double[x], i int[x], io int[x])';
[wlege, lused] = fmps(mex_id_, nlege, wlege, lw, lused, 1, lw, 1, 1);

