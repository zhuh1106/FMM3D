% mps test for clusters of sources: npts*n1^3
%
%

clear 
addpath('../../matlab/')

nd = 1; % num of charge vectors

n1 = 5;
nmpole = n1^3;
h = 1.0d0/(n1+1);
radius = h/8;
eps = 1e-15;

ntm = 16; % mterms, same mpole order for all mps centers. could change depending on complexity, closeness?
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
% set mps parameter
shift = h/10;
disp(['shift = ', num2str(shift), ' ']);
rscale = 1;
sc = shift;
if sc < 1 
  rscale = sc; % what is this doing?
end
% form all multipole expansion (just l3dformmpc, loop over cmpole)
[cmpole,rmpole,mterms,mpole,impole,ntot,nlege,wlege,lused] = myl3dformmpc(nd, rscale, source, charge, npts, ntm);
% multipole to local expansion
local = zeros(nd*ntot,1); ier = 0; 
local = lfmm3d_mps(nd, eps, nmpole, cmpole, rmpole, mterms, mpole, impole, local, ier);
% local expansion to pot
pot2 = myl3dtaevalp(nd, rmpole, cmpole, impole, local, mterms, source, npts, wlege, nlege);

disp(['mps vs direct error = ', num2str(max(abs(pot-pot2))/max(abs(pot))), ' ']);

figure(1),clf,
plot3(source(1,:),source(2,:),source(3,:),'.'); axis equal, hold on
plot3(cmpole(1,:),cmpole(2,:),cmpole(3,:),'*')


keyboard