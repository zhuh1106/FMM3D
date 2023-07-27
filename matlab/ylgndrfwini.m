function [w,lused] = ylgndrfwini(nmax, w, lw, lused)
% note: irat1=1, irat2=1+(nmax+1)^2 seem to be two starting indices 
% the fortran subroutine then calls ylgndrini, which initializes two size (nmax+1)^2 double precision arraries rat1 & rat2
% the two get assigned to w, which has a length of 2*(nmax+1)^2
% in the case of Laplace, l3dformmpc declares real *8 w(0:nlege,0:nlege), then calls ylgndrfw, which asks for w(irat1), w(irat2) when nterms (multipole expansion order) is greater than nlege
% Q: is rat2 for legendre derivative?
%
nw = length(w);
mex_id_ = 'ylgndrfwini(i int[x], io double[x], i int[x], io int[x])';
[w, lused] = fmm3d(mex_id_, nmax, w, lw, lused, 1, nw, 1, 1);
end

% ---------------------------------------------------------------------
