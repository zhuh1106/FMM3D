function [GSMat,shcMat,shvMat] = mykernelDLap(S, if_up, S_up)
% no type functionality...
% just return Galerkin SLP, did not change the function name, which is
% originally for Nystrom
%
if nargin < 2, if_up = false; end
if if_up % then assume S_up is provided, and p_up = p + 1;
  p = S.p;
  p_up = p+1; % upsampled nodes, assume p_up > p
  % build self eval mat (Nystrom on p_up grid)
  [SMat_up, SpMat, DMat] = kernelDLap(S_up); % Laplace SLP, SLPn, DLP self eval matrix
  shcMat_up = shAna(eye(2*p_up*(p_up+1))); % if even p, then abs(shcMat(end-2*p_up,:) - shcMat(end,:)) is 0, if odd, then + is 0 ...
  shvMat_up = shSyn(eye((p_up+1)^2),false);
  shcMat_up = shcMat_up(1:(p+1)^2,:); shvMat_up = shvMat_up(:,1:(p+1)^2); % val2coefsMat & coefs2valsMat (on p_up grid)
  L_up = interpsh(eye(2*p*(p+1)),p_up); % upsample interpolation matrix
  L_down = interpsh(eye(2*p_up*(p_up+1)),p); % downsample matrix
  GSMat = shcMat_up*(SMat_up)*shvMat_up; % Galerkin SLP mat on p
  shcMat = shcMat_up*L_up;
  shvMat = L_down*shvMat_up;

else
  [SMat, SpMat, DMat] = kernelDLap0(S);

  keyboard
end

end