% not quite mps...
% 
% 1. gmres call function that maps spherical harmonics coeffcients to spherical harmonics coefficients via mps
% 2. need to extend to ~100 particles...
%
% Hai 07/30/23, finalize laplace mps interface 07/31/23

profile clear
profile on
clear all

addpath ./spharm/; 
addpath ../../matlab/;

p = 16; % 3d spherical harmonics order
[S0,S0_up] = mySurfaceSph(p);


%
%===== build self eval mat, same for all =====
if_up = true;
[GSMat0,shcMat,shvMat] = mykernelDLap(S0, if_up, S0_up); % Galerkin SLP mat, name is misleading (follow Nystrom) 
                                                         % now full rank, temp fix, user unaware


% 
%===== follow mps test geometry configuration =====
n1 = 3; % ~ 30s, 38 iteration, 8 core macbook pro
% n1 = 5; % ~ 3min, 53 iteration, 8 core macbook pro
h = 1.0d0/(n1+1);
radius = h/8;
s = []; y_source.x = []; y_source.w = []; pt_source = [];
ijk = 0;
for i = 1:n1
  for j = 1:n1 
    for k = 1:n1 
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
end


% 
%===== boundary condition =====
sx = cellfun(@(p)p.x,s,'uniformoutput',0); snx = cellfun(@(p)p.nx,s,'uniformoutput',0); sw = cellfun(@(p)p.w,s,'uniformoutput',0); 
sx = horzcat(sx{:}); snx = horzcat(snx{:}); sw = horzcat(sw{:});
S = struct('x',sx,'nx',snx,'w',sw);
rhs = Lap3dSLPfmm(S,y_source,pt_source(:),1e-14);
rhsc = zeros(numel(s)*(p+1)^2,1);
for k=1:numel(s) % from rhs to rhs coefs
  rhsc((k-1)*(p+1)^2+(1:(p+1)^2)) = shcMat*rhs((k-1)*((p+1)*2*p)+(1:(p+1)*2*p));
end


% 
%===== mps setup =====
% for laplace mps call
nd = 1; rscale = radius; 
npts = (p+1)*2*p; ntm = p;
source = S.x; nmpole = numel(s);
cmpole = zeros(3,nmpole);
rmpole = zeros(nmpole,1); % rescaling factor for each mps center
mterms = zeros(nmpole,1); % multipole expansion order 
impole = zeros(nmpole,1); % index of mpole for each mps center
ntot = 0; % total num of multipole expansion coefficients
impole(1) = 1;
for i = 1:nmpole
  rmpole(i) = rscale;
  mterms(i) = ntm;
  cmpole(1,i) = mean(source(1,((i-1)*npts+1):i*npts)); % maybe should choose average as mpole center
  cmpole(2,i) = mean(source(2,((i-1)*npts+1):i*npts));
  cmpole(3,i) = mean(source(3,((i-1)*npts+1):i*npts));
  ilen = (mterms(i)+1)*(2*mterms(i)+1);
  ntot = ntot + (mterms(i)+1)*(2*mterms(i)+1);
  if i<nmpole, impole(i+1) = impole(i) + nd*ilen; end
end
% for mpole/local & sh coefs map
[GSMatn,local_shc1n,local_shc2n,mpole_shv1n,mpole_shv2n] = lfmm3d_mps_setup(GSMat0,shcMat,shvMat,s,radius,p);


% 
%===== mps solve =====
myoperator = @(x) lfmm3d_mps_op(GSMatn,local_shc1n,local_shc2n,mpole_shv1n,mpole_shv2n,x,nd,nmpole,cmpole,rmpole,mterms,impole);
mucS = gmres(@(x) myoperator(x),rhsc,[],1e-14,100); % no null space
muS = zeros(numel(s)*(p+1)*2*p,1);
for k=1:numel(s) % from slp mu coefs to mu
  muS((k-1)*((p+1)*2*p)+(1:(p+1)*2*p)) = real(shvMat*mucS((k-1)*(p+1)^2+(1:(p+1)^2)));
end

% 
%===== verify solution =====
t.x = [S.x + radius*S.nx, S.x + radius*3/4*S.nx, S.x + radius*1/2*S.nx];
fhom = Lap3dSLPfmm(t,y_source,pt_source(:),1e-14);
fnumS = Lap3dSLPfmm(t,S,muS(:),1e-14);
errS = abs(fnumS-fhom)/max(abs(fhom));
% figure(1),clf,scatter3(S.x(1,:),S.x(2,:),S.x(3,:),[],rhs,'.'); hold on, axis equal, plot3(t.x(1,:),t.x(2,:),t.x(3,:),'.')
figure(2),clf,
scatter3(t.x(1,:),t.x(2,:),t.x(3,:),10,log10(errS),'filled');
colorbar, axis equal, hold on
for j=1:numel(s)
  surf(reshape(s{j}.x(1,:),p+1,[]),reshape(s{j}.x(2,:),p+1,[]),reshape(s{j}.x(3,:),p+1,[]))
end
caxis([-12 0])
title(['Galerkin SLP p=' num2str(p) ', ' num2str(numel(rhsc)) 'x' num2str(numel(rhsc)) ', rank ' num2str(rank(GSMat0)*numel(s)) ', max err ' num2str(max(errS))], 'FontSize', 16)

profile viewer

keyboard