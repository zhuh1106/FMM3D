function [evecs,hvecs]=em3taevaleh(rk,center,aimpole,bimpole,nterms, ...
    sphere_xyz,sphere_r, ...
    rnodes,weights,nphi,ntheta)
%em3taevaleh:  Evaluate incoming EM multipole on a spherical grid.
%
%  [evecs,hvecs]=em3taevaleh(rk,center,aimpole,bimpole,nterms,...
%         sphere_xyz,sphere_r,rnodes,weights,nphi,ntheta)
%
%  All EM multipoles are (NTERMS+1)-by-(2*NTERMS+1) complex matrices.
%
%  Input parameters:
%
%    rk - the frequency parameter
%    center - center of the multipole expansion
%    aimpole,bimpole - complex(ncoefs): incoming EM multipoles
%    nterms - the number of terms in multipole expansions 
%    sphere_xyz - real(3): sphere center location (for E and H fields)
%    sphere_r - sphere radius (for E and H fields)
%    rnodes,weights,nphi,ntheta - 
%          spherical grid, constructed via a preceding call to e3fgrid(nterms)
%                 
%  Output parameters:
%
%    evecs,hvecs - E and H values on the spherical grid.
%
%

nnodes=nphi*ntheta;
evecs = zeros(3,nnodes) + 1i*zeros(3,nnodes);
hvecs = zeros(3,nnodes) + 1i*zeros(3,nnodes);

mex_id_ = 'em3taevaleh(i dcomplex[x], i double[], i dcomplex[], i dcomplex[], i int[x], i double[], i double[], io dcomplex[], io dcomplex[], i double[], i double[], i int[x], i int[x])';
[evecs, hvecs] = fmps(mex_id_, rk, center, aimpole, bimpole, nterms, sphere_xyz, sphere_r, evecs, hvecs, rnodes, weights, nphi, ntheta, 1, 1, 1, 1);


