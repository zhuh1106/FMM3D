% 1. fix the null space issue via upsampling... maybe p -> p+1 or p+2
% then compute Nystrom, then downsample to check rank
%
% 2. solve a multi particle bvp in coefficient space
%
% Hai 07/27/23

clear all

addpath ./spharm/; 

p = 24; % 3d spherical harmonics order
[S0,S0_up] = mySurfaceSph(p);

% build self eval mat, same for all
if_up = true;
[GSMat0,shcMat,shvMat] = mykernelDLap(S0, if_up, S0_up); % Galerkin SLP mat, name is misleading (follow Nystrom)
                                                      % now full rank, temp fix, user unaware

% follow mps test geometry configuration
n1 = 2;
h = 1.0d0/(n1+1);
radius = h/8;
s = []; y_source.x = []; y_source.w = []; pt_source = [];
ijk = 0;
for i = 1:n1
  for j = 1:n1 
    k = 1; % add loop later
      ijk = ijk + 1;
      s{ijk} = S0; % copy
      s{ijk}.x(1,:) = h*i + radius*S0.x(1,:); % scale and shift x
      s{ijk}.x(2,:) = h*j + radius*S0.x(2,:);
      s{ijk}.x(3,:) = h*k + radius*S0.x(3,:);
      s{ijk}.w = radius^2*S0.w; % scale quadrature weights
      % for boundary condition setup
      y_source.x = [y_source.x h*[i;j;k]+radius/2*(rand(3,1)-1/2)];
      y_source.w = [y_source.w 1];
      pt_source = [pt_source 1];
  end
end

% build off diagonal block (hard-coded only on 2x2x1 for now)
GSMat0 = radius*GSMat0; % rescale SLP diagonal block
GSMat12 = shcMat*Lap3dSLPmat(s{1},s{2})*shvMat; GSMat13 = shcMat*Lap3dSLPmat(s{1},s{3})*shvMat; GSMat14 = shcMat*Lap3dSLPmat(s{1},s{4})*shvMat;
GSMat21 = shcMat*Lap3dSLPmat(s{2},s{1})*shvMat; GSMat23 = shcMat*Lap3dSLPmat(s{2},s{3})*shvMat; GSMat24 = shcMat*Lap3dSLPmat(s{2},s{4})*shvMat;
GSMat31 = shcMat*Lap3dSLPmat(s{3},s{1})*shvMat; GSMat32 = shcMat*Lap3dSLPmat(s{3},s{2})*shvMat; GSMat34 = shcMat*Lap3dSLPmat(s{3},s{4})*shvMat;
GSMat41 = shcMat*Lap3dSLPmat(s{4},s{1})*shvMat; GSMat42 = shcMat*Lap3dSLPmat(s{4},s{2})*shvMat; GSMat43 = shcMat*Lap3dSLPmat(s{4},s{3})*shvMat;

% build entire system matrix coeffs2coeffs
GSMat = [ GSMat0   GSMat12   GSMat13   GSMat14;...
         GSMat21    GSMat0   GSMat23   GSMat24;...
         GSMat31   GSMat32    GSMat0   GSMat34;...
         GSMat41   GSMat42   GSMat43    GSMat0];

% boundary condition
S = struct('x',[s{1}.x s{2}.x s{3}.x s{4}.x],'nx',[s{1}.nx s{2}.nx s{3}.nx s{4}.nx],'w',[s{1}.w s{2}.w s{3}.w s{4}.w]);
Ahom = Lap3dSLPmat(S,y_source);
rhs = Ahom*pt_source(:); 
rhsc = [shcMat*rhs(1:end/4);shcMat*rhs(end/4+1:2*end/4);shcMat*rhs(2*end/4+1:3*end/4);shcMat*rhs(3*end/4+1:end)];

% solve
mucS = GSMat\(rhsc(:)); % no null space
muS = real([shvMat*mucS(1:end/4);shvMat*mucS(end/4+1:2*end/4);shvMat*mucS(2*end/4+1:3*end/4);shvMat*mucS(3*end/4+1:end)]); % coefs2vals
% figure(1),plot(log10(abs(muS)),'.')

% verify solution
t.x = [S.x + radius*S.nx, S.x + radius*3/4*S.nx, S.x + radius*1/2*S.nx];
fhom = Lap3dSLPmat(t,y_source)*pt_source(:);
fnumS = Lap3dSLPmat(t,S)*muS;
errS = abs(fnumS-fhom)/max(abs(fhom));
% figure(1),clf,scatter3(S.x(1,:),S.x(2,:),S.x(3,:),[],rhs,'.'); hold on, axis equal, plot3(t.x(1,:),t.x(2,:),t.x(3,:),'.')

figure(2),clf,
scatter3(t.x(1,:),t.x(2,:),t.x(3,:),10,log10(errS),'filled');
colorbar, axis equal, hold on
for j=1:numel(s)
  surf(reshape(s{j}.x(1,:),p+1,[]),reshape(s{j}.x(2,:),p+1,[]),reshape(s{j}.x(3,:),p+1,[]))
end
caxis([-12 0])
title(['Galerkin SLP p=' num2str(p) ', ' num2str(size(GSMat,1)) 'x' num2str(size(GSMat,2)) ', rank ' num2str(rank(GSMat)) ', max err ' num2str(max(errS))], 'FontSize', 16)

%%% let's verify a few things before a general solve
%%% assume we want to verify lhs = rhs in spherical harmonics expansion
% 1st check 1st row, without diagonal block, SLPmat vs SLPfmm (check)
GSMat1 = [ zeros(size(GSMat0))   GSMat12   GSMat13   GSMat14]; % 1st row
SMat1 = [ zeros(size(Lap3dSLPmat(s{1},s{1}))) Lap3dSLPmat(s{1},s{2}) Lap3dSLPmat(s{1},s{3}) Lap3dSLPmat(s{1},s{4})];
u1 = SMat1*muS;
S1 = struct('x',[s{2}.x s{3}.x s{4}.x],'nx',[s{2}.nx s{3}.nx s{4}.nx],'w',[s{2}.w s{3}.w s{4}.w]);
u1_fmm = Lap3dSLPfmm(s{1},S1,muS(end/4+1:end),eps);
max(abs(u1-u1_fmm))
% 2nd swtich to a laplace mps call
nd = 1; rscale = radius; % is this a good rscale
source = [s{1}.x s{2}.x s{3}.x s{4}.x]; charge = muS'.*[s{1}.w s{2}.w s{3}.w s{4}.w];
npts = (p+1)*2*p; ntm = p; nmpole = numel(s);
[cmpole,rmpole,mterms,mpole,impole,ntot,nlege,wlege,lused] = myl3dformmpc(nd, rscale, source, charge, npts, ntm);
local = zeros(nd*ntot,1); ier = 0; 
local = lfmm3d_mps(nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);
pot = myl3dtaevalp(nd, rmpole, cmpole, impole, local, mterms, source, npts, wlege, nlege);
u1_mps = 1/(4*pi)*pot(1:end/4)';
max(abs(u1-u1_mps))
% 3rd combine diagonal block & mps off diagonal
shcMatcell = repmat({shcMat}, 1, numel(s));
u_mps = 1/(4*pi)*pot(:);
uhat_mps = blkdiag(shcMatcell{:})*u_mps;
GSMatcell = repmat({GSMat0}, 1, numel(s));
uhat_diag = blkdiag(GSMatcell{:})*mucS;
uhat = uhat_mps + uhat_diag;

keyboard