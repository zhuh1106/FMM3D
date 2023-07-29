% not quite mps...
% 
% 1. gmres call function that maps spherical harmonics coeffcients to spherical harmonics coefficients via mps
%
% 2. figure out a way to switch to a different set of coefficients (need to be careful with dimension, rank for different coeffs sets...)
%
% Hai 07/28/23

clear all

addpath ./spharm/; 
addpath ../../matlab/;

p = 16; % 3d spherical harmonics order
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
if 0 % form dense operator
GSMat_operator = @(x) mylfmm3d_mps(GSMat,x);
mucS = gmres(@(x) GSMat_operator(x),rhsc,[],1e-14,100); % no null space
end
shcMatcell = repmat({sparse(shcMat)}, 1, numel(s)); shcMatall = blkdiag(shcMatcell{:});
shvMatcell = repmat({sparse(shvMat)}, 1, numel(s)); shvMatall = blkdiag(shvMatcell{:});
GSMatcell = repmat({sparse(GSMat0)}, 1, numel(s)); GSMatall = blkdiag(GSMatcell{:});
nd = 1; rscale = radius; s; p;
if 0
GSMat_operator2 = @(x) mylfmm3d_mps2(GSMat0,shcMatall,shvMatall,x,nd,rscale,s,p);
mucS = gmres(@(x) GSMat_operator2(x),rhsc,[],1e-14,100); % no null space
end
GSMat_operator3 = @(x) mylfmm3d_mps3(GSMat0,shcMatall,shvMatall,x,nd,rscale,s,p);
mucS = gmres(@(x) GSMat_operator3(x),rhsc,[],1e-14,100); % no null space
muS = real([shvMat*mucS(1:end/4);shvMat*mucS(end/4+1:2*end/4);shvMat*mucS(2*end/4+1:3*end/4);shvMat*mucS(3*end/4+1:end)]); % coefs2vals

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
title(['Galerkin SLP p=' num2str(p) ', ' num2str(numel(rhsc)) 'x' num2str(numel(rhsc)) ', rank ' num2str(rank(GSMat)) ', max err ' num2str(max(errS))], 'FontSize', 16)


keyboard

function u = mylfmm3d_mps(GSMat,rhs)
u = GSMat*rhs;
end

function u = mylfmm3d_mps2(GSMat0,shcMat,shvMat,rhshat,nd,rscale,s,p)
% imaginary part does matter!!! 07/28/23
% matrix version
if 0
SMat12 = Lap3dSLPmat(s{1},s{2}); SMat13 = Lap3dSLPmat(s{1},s{3}); SMat14 = Lap3dSLPmat(s{1},s{4});
SMat21 = Lap3dSLPmat(s{2},s{1}); SMat23 = Lap3dSLPmat(s{2},s{3}); SMat24 = Lap3dSLPmat(s{2},s{4});
SMat31 = Lap3dSLPmat(s{3},s{1}); SMat32 = Lap3dSLPmat(s{3},s{2}); SMat34 = Lap3dSLPmat(s{3},s{4});
SMat41 = Lap3dSLPmat(s{4},s{1}); SMat42 = Lap3dSLPmat(s{4},s{2}); SMat43 = Lap3dSLPmat(s{4},s{3});

SMat = [ zeros(size(SMat12))   SMat12   SMat13   SMat14;...
         SMat21    zeros(size(SMat12))  SMat23   SMat24;...
         SMat31   SMat32    zeros(size(SMat12))  SMat34;...
         SMat41   SMat42   SMat43    zeros(size(SMat12))];
u_mps = shcMat*SMat*(shvMat*rhshat);
uhat_diag = blkdiag(GSMat0,GSMat0,GSMat0,GSMat0)*rhshat(:);
u = uhat_diag + u_mps;
end
% keyboard
% fmm version 
rhsr = real(shvMat*rhshat); rhsi = imag(shvMat*rhshat);
S1 = struct('x',[s{2}.x s{3}.x s{4}.x],'nx',[s{2}.nx s{3}.nx s{4}.nx],'w',[s{2}.w s{3}.w s{4}.w]);
u1_fmmr = Lap3dSLPfmm(s{1},S1,rhsr(end/4+1:end),eps); % 1st row
u1_fmmi = Lap3dSLPfmm(s{1},S1,rhsi(end/4+1:end),eps);
S2 = struct('x',[s{1}.x s{3}.x s{4}.x],'nx',[s{1}.nx s{3}.nx s{4}.nx],'w',[s{1}.w s{3}.w s{4}.w]);
u2_fmmr = Lap3dSLPfmm(s{2},S2,[rhsr(1:end/4);rhsr(end/2+1:end)],eps); % 2nd row
u2_fmmi = Lap3dSLPfmm(s{2},S2,[rhsi(1:end/4);rhsi(end/2+1:end)],eps);
S3 = struct('x',[s{1}.x s{2}.x s{4}.x],'nx',[s{1}.nx s{2}.nx s{4}.nx],'w',[s{1}.w s{2}.w s{4}.w]);
u3_fmmr = Lap3dSLPfmm(s{3},S3,[rhsr(1:end/2);rhsr(3*end/4+1:end)],eps); % 3rd row
u3_fmmi = Lap3dSLPfmm(s{3},S3,[rhsi(1:end/2);rhsi(3*end/4+1:end)],eps);
S4 = struct('x',[s{1}.x s{2}.x s{3}.x],'nx',[s{1}.nx s{2}.nx s{3}.nx],'w',[s{1}.w s{2}.w s{3}.w]);
u4_fmmr = Lap3dSLPfmm(s{4},S4,rhsr(1:3*end/4),eps); % 3rd row
u4_fmmi = Lap3dSLPfmm(s{4},S4,rhsi(1:3*end/4),eps);
uhat_diag = blkdiag(GSMat0,GSMat0,GSMat0,GSMat0)*rhshat(:);
u = uhat_diag + shcMat*[u1_fmmr+1i*u1_fmmi;u2_fmmr+1i*u2_fmmi;u3_fmmr+1i*u3_fmmi;u4_fmmr+1i*u4_fmmi];

end

function uhat = mylfmm3d_mps3(GSMat,shcMat,shvMat,rhshat,nd,rscale,s,p)
% demo case

%%% mps 
rhsr = real(shvMat*rhshat); rhsi = imag(shvMat*rhshat);
source = [s{1}.x s{2}.x s{3}.x s{4}.x]; 
npts = (p+1)*2*p; ntm = p; nmpole = numel(s); ier = 0; 
% real
charger = rhsr.'.*[s{1}.w s{2}.w s{3}.w s{4}.w];
[cmpole,rmpole,mterms,mpole,impole,ntot,nlege,wlege,lused] = myl3dformmpc(nd, rscale, source, charger, npts, ntm);
local = zeros(nd*ntot,1); 
local = lfmm3d_mps(nd, 1e-16, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);
potr = myl3dtaevalp(nd, rmpole, cmpole, impole, local, mterms, source, npts, wlege, nlege);
uhat_mpsr = 1/(4*pi)*shcMat*potr(:);
% imaginary
chargei = rhsi.'.*[s{1}.w s{2}.w s{3}.w s{4}.w];
[cmpole,rmpole,mterms,mpole,impole,ntot,nlege,wlege,lused] = myl3dformmpc(nd, rscale, source, chargei, npts, ntm);
local = zeros(nd*ntot,1); 
local = lfmm3d_mps(nd, 1e-16, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);
poti = myl3dtaevalp(nd, rmpole, cmpole, impole, local, mterms, source, npts, wlege, nlege);
uhat_mpsi = 1/(4*pi)*shcMat*poti(:);
% uhat
uhat_mps = uhat_mpsr + 1i*uhat_mpsi;

%%% self
uhat_diag = blkdiag(GSMat,GSMat,GSMat,GSMat)*rhshat(:);

%%% together
uhat = uhat_diag + uhat_mps;

end


