function uhat = lfmm3d_mps_op(GSMatn,local_shc1n,local_shc2n,mpole_shv1n,mpole_shv2n,rhshat,nd,nmpole,cmpole,rmpole,mterms,impole)
% laplace mps operator fun
% more or less finalize the interface
% this one finishes sh coefs to mpole map & local to sh coefs map...
% lfmm3d_mps calls fortran mex
%
% gmres should iterate on this operator
% within this operator
% 1st  'rhshat' (sh coefs) --> 'mpole' (multipole expansion)
% then 'mpole' --> 'local' (Libin's fortran Laplace mps call)
% then 'local' (local expansion) --> 'uhat' (sh coefs)
%
% diagonal part is just 'rhshat' (sh coefs) --> 'uhat' sh coefs 
% (maybe work on a idensity diagonal block version.. 
%  similar to Zydrunas' fmps?)
% (or get rid of the middle man, a proper mps interface.. 
%  similar to Jun, Motoki and Leslie's A fast solver for multi-particle scattering in a layered medium)
%
% extend to Stokes...?
%
% Hai 07/31/23

%===== off-diagonal part (mps) =====
ntot = sum((mterms+1).*(2*mterms+1)); ier = 0; % ntm = p % parameter setup to call lfmm3d_mps
%----- real part of rhshat --------
% sh coefs to mpole coefs
mpoler = mpole_shv1n*real(rhshat) + mpole_shv2n*imag(rhshat);
% m2l
localr = zeros(nd*ntot,1); % maybe use nd = 2 later... or combine real & imaginary part
localr = lfmm3d_mps(nd, 1e-16, nmpole, cmpole, rmpole, mterms, mpoler, impole, localr, ier); % call fortran lfmm3d_mps
% l2pot
uhat_mpsr = 1/(4*pi)*(local_shc1n*real(localr)+local_shc2n*imag(localr));
%----- imaginary part of rhshat --------
% sh coefs to mpole coefs
mpolei = -mpole_shv2n*real(rhshat) + mpole_shv1n*imag(rhshat);
% m2l
locali = zeros(nd*ntot,1); 
locali = lfmm3d_mps(nd, 1e-16, nmpole, cmpole, rmpole, mterms, mpolei, impole, locali, ier); % call fortran lfmm3d_mps
% l2pot
uhat_mpsi = 1/(4*pi)*(local_shc1n*real(locali)+local_shc2n*imag(locali));
%----- uhat = real + imaginary --------
uhat_mps = uhat_mpsr + 1i*uhat_mpsi;

%===== diagonal part (sh coefs) =====
uhat_diag = GSMatn*rhshat(:);

%===== 1 mat-vec product in sh coefs space =====
uhat = uhat_diag + uhat_mps;

end