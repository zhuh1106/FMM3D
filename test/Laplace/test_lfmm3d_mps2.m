% mps test for clusters of sources: npts*n1^3
%
%

profile clear
profile on

clear 
addpath('../../matlab/')

nd = 1; % num of charge vectors

n1 = 5;
nmpole = n1^3;
h = 1.0d0/(n1+1);
radius = h/5;

nt = 19;
ntm = 8; % mterms, same mpole order for all mps centers. could change depending on complexity, closeness?
          % if change source to ~(rand(1)-1/2)*2*h/10, need to increase ntm to ~32
S = mySurfaceSph(ntm); % 3d spherical harmonics
npts = numel(S.x(1,:)); % num of sources per particle, assume same spherical harmonics discretization
ns = npts*n1^3;

% generate sources uniformly in the unit cube 
source=zeros(3,ns);
charge=zeros(nd,ns);
ijk = 0;
for i = 1:n1
  for j = 1:n1
    for k = 1:n1
      for ell = 1:npts
        ijk = ijk + 1;
        source(1,ijk) = h*i + radius*S.x(1,ell);
        source(2,ijk) = h*j + radius*S.x(2,ell);
        source(3,ijk) = h*k + radius*S.x(3,ell);
      end
    end
  end
end
dnorm = 0;
for i=1:ns
  for idim=1:nd
    charge(idim,i) = i/ns;
    charge(idim,i) = hkrand_mex(0);
    dnorm = dnorm + abs(charge(idim,i))^2;
  end
end
dnorm = sqrt(dnorm);
for i=1:ns
  for idim = 1:nd
    charge(idim,i) = charge(idim,i)/dnorm;
  end
end 

% do the direct calculation (source to source)
thresh = 1.0d-15;
ifcharge = 1; ifdipole = 0; ifpgh = 1;
ifpghtarg = 0; iper = 0;
srcinfo.sources = source;
srcinfo.charges = charge;
srcinfo.nd = nd;
U1 = lfmm3d(thresh,srcinfo,1);
pot = U1.pot;
for i = 1:nmpole % loop over mps center, subtract particle self
  s = []; s.x = srcinfo.sources(:,((i-1)*npts+1):i*npts); s.w = ones(1,npts);
  A = Lap3dSLPmat(s,s); A(1:npts+1:end) = 0;
  pot(((i-1)*npts+1):i*npts) = pot(((i-1)*npts+1):i*npts) - 4*pi*srcinfo.charges(((i-1)*npts+1):i*npts)*A';
end

%%% mps setup
% set mps center (center of each cluster of particles)
cmpole = zeros(3,nmpole);
shift = h/10;
disp(['shift = ', num2str(shift), ' ']);
for i = 1:nmpole
  cmpole(1,i) = mean(source(1,((i-1)*npts+1):i*npts)); % maybe should choose average as mpole center
  cmpole(2,i) = mean(source(2,((i-1)*npts+1):i*npts));
  cmpole(3,i) = mean(source(3,((i-1)*npts+1):i*npts));
end
% form a multipole expansion at each center
mterms=zeros(nmpole,1);
impole=zeros(nmpole,1); % index of mpole for each mps center
ntot = 0; % total num of multipole expansion coefficients
for i = 1:nmpole
  mterms(i) = ntm;
  ntot = ntot + (mterms(i)+1)*(2*mterms(i)+1);
end
mpole = zeros(nd*ntot,1); % nd = 1 for now, complex
impole(1) = 1;
for i = 1:nmpole-1 % loop over mps centers
  ilen = (mterms(i)+1)*(2*mterms(i)+1);
  impole(i+1) = impole(i) + nd*ilen;
end 
nlege = 300;
lw = 2*(nlege+1)^2; % dimension
wlege = zeros(lw,1);
lused = 0;
[wlege,lused] = ylgndrfwini_mex(nlege, wlege, lw, lused); % recursion coefficients
npts; % 
rscale = 1;
sc = shift;
if sc < 1 
  rscale = sc; % what is this doing?
end
rmpole = zeros(nmpole,1); % rescaling factor for each mps center
for i = 1:nmpole
  rmpole(i) = rscale;
  mpoletmp = zeros(nd,mterms(i)+1,2*mterms(i)+1);
  mpoletmp =  l3dformmpc_mex(nd, rscale, source(:,((i-1)*npts+1):i*npts), charge(:,((i-1)*npts+1):i*npts), npts, cmpole(:,i), mterms(i), mpoletmp, wlege, nlege);
  mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) = mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) + mpoletmp(:);
end % why imaginary part so small?
% mps call: multipole to local operator
local = zeros(nd*ntot,1); ier = 0; 
local = lfmm3d_mps0_mex(nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);

% post process
interms = (mterms(1)+1)*(2*mterms(1)+1); % same for all mps centers
pot2 = zeros(nd,ns);
npts;
for i = 1:nmpole % loop over mps center, ship local expansion to targ
  pot2tmp = zeros(nd,npts);
  pot2tmp = l3dtaevalp_mex(nd, rmpole(i), cmpole(:,i), local(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1), mterms(i), source(:,((i-1)*npts+1):i*npts), npts, pot2tmp, wlege, nlege);
  pot2(:,((i-1)*npts+1):i*npts) = pot2(:,((i-1)*npts+1):i*npts) + pot2tmp;
end

disp(['mps vs direct error = ', num2str(max(abs(pot-pot2))/max(abs(pot))), ' ']);

figure(1),clf,
plot3(source(1,:),source(2,:),source(3,:),'.'); axis equal, hold on
plot3(cmpole(1,:),cmpole(2,:),cmpole(3,:),'*')

profile viewer

keyboard