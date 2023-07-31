function [GSMatn,local_shc1n,local_shc2n,mpole_shv1n,mpole_shv2n] = lfmm3d_mps_setup(GSMat0,shcMat,shvMat,s,radius,p)
% additional input required for maps between mpole/local and middle man (sh coefs)
% other parameter setup done outside
% Stokes maybe also be doable this way... need to take a closer look later
%
% GSMat0 is for unit sphere, need to rescale...
% 
% Hai 07/31/23

nd = 1; rscale = radius; 
% here build map from sh coefs to mpole on each particle
rhs_mat = eye((p+1)^2);
mpole_shvtmp1 = zeros((p+1)*(2*p+1),(p+1)^2); % for real & imaginary parts of rhshat (mpole <-> shv)
mpole_shvtmp2 = zeros((p+1)*(2*p+1),(p+1)^2); % 1 or 2 does not mean anything, maybe come up with a better name
npts = (p+1)*2*p; ntm = p;
for k=1:(p+1)^2
  [~,~,~,mpole1,~,~,~,~,~] = myl3dformmpc(nd, rscale, s{1}.x, real(shvMat*rhs_mat(:,k)).'.*s{1}.w, npts, ntm);
  [~,~,~,mpole2,~,~,~,~,~] = myl3dformmpc(nd, rscale, s{1}.x, real(shvMat*1i*rhs_mat(:,k)).'.*s{1}.w, npts, ntm);
  mpole_shvtmp1(:,k) = mpole1;
  mpole_shvtmp2(:,k) = mpole2;
end
mpole_shv1 = sparse((p+1)*(2*p+1),(p+1)^2); mpole_shv1(abs(mpole_shvtmp1)>1e-14) = mpole_shvtmp1(abs(mpole_shvtmp1)>1e-14); % this sparse pattern should be derived on paper ...
mpole_shv2 = sparse((p+1)*(2*p+1),(p+1)^2); mpole_shv2(abs(mpole_shvtmp2)>1e-14) = mpole_shvtmp2(abs(mpole_shvtmp2)>1e-14);
mpole_shv1_cell = repmat({sparse(mpole_shv1)}, 1, numel(s)); mpole_shv2_cell = repmat({sparse(mpole_shv2)}, 1, numel(s)); 
mpole_shv1n = blkdiag(mpole_shv1_cell{:}); mpole_shv2n = blkdiag(mpole_shv2_cell{:}); % blkdiag matrix that repeats single particle sh coefs to mpole map
% here build map from local to sh coefs on each particle
idx =reshape(1:(2*p+1)*(p+1), p+1, 2*p+1);
for ii=0:p, idx(ii+1, [1:p-ii p+ii+2:2*p+1]) = 0; end
idx = idx'; idx = idx(idx~=0); 
idxmat = zeros((p+1)^2,(2*p+1)*(p+1)); % to expand triangular local expansion to rectangular expansion
for ii=1:(p+1)^2, idxmat(ii,idx(ii)) = 1; end
localmat = eye((p+1)^2);
potmat1 = zeros(2*p*(p+1),(p+1)^2);
potmat2 = zeros(2*p*(p+1),(p+1)^2);
nlege = 300; lw = 2*(nlege+1)^2; wlege = zeros(lw,1); lused = 0;
[wlege,lused] = ylgndrfwini(nlege, wlege, lw, lused); % recursion coefficients
for k=1:(p+1)^2
  localk = localmat(:,k); 
  localexpandk = idxmat'*localk;
  potk = myl3dtaevalp(nd, rscale, mean(s{1}.x,2), 1, localexpandk.', p, s{1}.x, npts, wlege, nlege);
  potmat1(:,k) = potk;
  localexpandk = idxmat'*(1i*localk); 
  potk = myl3dtaevalp(nd, rscale, mean(s{1}.x,2), 1, localexpandk.', p, s{1}.x, npts, wlege, nlege);
  potmat2(:,k) = potk;
end
local_shctmp1 = shcMat*potmat1*idxmat; local_shctmp2 = shcMat*potmat2*idxmat; % dense ones
local_shc1 = sparse((p+1)^2,(p+1)*(2*p+1)); local_shc1(abs(local_shctmp1)>1e-14) = local_shctmp1(abs(local_shctmp1)>1e-14); % this sparse pattern should be derived on paper ...
local_shc2 = sparse((p+1)^2,(p+1)*(2*p+1)); local_shc2(abs(local_shctmp2)>1e-14) = local_shctmp2(abs(local_shctmp2)>1e-14); 
local_shc1_cell = repmat({sparse(local_shc1)}, 1, numel(s)); local_shc2_cell = repmat({sparse(local_shc2)}, 1, numel(s)); 
local_shc1n = blkdiag(local_shc1_cell{:}); local_shc2n = blkdiag(local_shc2_cell{:}); % blkdiag matrix that repeats single particle local to sh coefs map
%
GSMat0 = radius*GSMat0; GSMatcell = repmat({sparse(GSMat0)}, 1, numel(s)); GSMatn = blkdiag(GSMatcell{:});

end