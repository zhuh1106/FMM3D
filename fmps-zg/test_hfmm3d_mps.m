%
% run make -f makefile_macosx in this folder
% run make -f test_hfmm3d_mps.make in /test/Helmholtz folder

done = 1;
pi = 4*atan(done);

nd = 3;


n1 = 3;
ns = n1^3;
nc = ns;

nt = 19;
eps = 0.5e-9;

boxsize = 1;
dlam = 1/boxsize;
zk = 2*pi/dlam + 1i*0.02;

% generate sources uniformly in the unit cube 
h = 1.0d0/(n1+1);
ijk = 0;
source = zeros(3,ns);
for i = 1:n1
  for j = 1:n1
    for k = 1:n1
      ijk = ijk + 1;
      source(1,ijk) = h*i;
      source(2,ijk) = h*j;
      source(3,ijk) = h*k;
    end
  end
end

dnorm = 0;
charge = zeros(nd,ns);
dipvec = zeros(nd,3,ns);
pot = zeros(nd,ns);
grad = zeros(nd,3,ns);
for i=1:ns
  for idim=1:nd
    charge(idim,i) = hkrand2(0) + 1i*hkrand2(0);
    dnorm = dnorm + abs(charge(idim,i))^2;

    dipvec(idim,1,i) = hkrand2(0) + 1i*hkrand2(0);
    dipvec(idim,2,i) = hkrand2(0) + 1i*hkrand2(0);
    dipvec(idim,3,i) = hkrand2(0) + 1i*hkrand2(0);
  end
end

dnorm = sqrt(dnorm);
for i=1:ns
  for idim = 1:nd
    charge(idim,i) = charge(idim,i)/dnorm;
  end
end
charge(1,:) = [0.16710E-01  0.48855E-01  0.96766E-01  0.61500E-01  0.11416E+00  0.42110E-01 ...
  0.11048E-01  0.51491E-01  0.74365E-01  0.15395E-01  0.82934E-01  0.21615E-01 ...
  0.38634E-01  0.24700E-01  0.10281E+00  0.11201E+00  0.12970E+00  0.12977E+00 ...
  0.19869E-01  0.99873E-01  0.38627E-01  0.87639E-01  0.10591E+00  0.12004E+00 ...
  0.11179E+00  0.61330E-01  0.12556E+00] ...
  +1i *[0.89192E-01  0.15050E-01  0.65964E-01  0.23289E-01  0.54856E-01  0.76195E-01 ...
  0.41835E-01  0.97569E-01  0.74882E-01  0.38757E-01  0.29141E-01  0.93632E-01 ...
  0.11209E+00  0.10488E+00  0.15659E-01  0.16592E-01  0.24828E-01  0.42027E-01 ...
  0.50221E-01  0.46430E-01  0.32764E-01  0.18467E-04  0.42915E-01  0.31371E-01 ...
  0.59298E-01  0.13532E+00  0.10426E+00]; 
charge(2,:) = [0.10632E+00  0.10227E+00  0.94702E-01  0.17474E-01  0.92729E-01  0.15694E-02 ...
  0.61361E-01  0.12888E+00  0.29222E-01  0.40251E-01  0.43134E-01  0.13391E+00 ...
  0.65239E-01  0.76331E-01  0.71831E-03  0.11844E-01  0.40062E-01  0.14298E-01 ...
  0.11460E+00  0.11195E+00  0.34582E-01  0.30574E-01  0.78154E-01  0.16348E-01 ...
  0.13735E+00  0.82523E-01  0.98532E-01] ...
  +1i *[0.11731E+00  0.12237E+00  0.64095E-02  0.55654E-02  0.12165E-01  0.79841E-01 ...
  0.96316E-01  0.75861E-01  0.78068E-01  0.71808E-01  0.75239E-02  0.45521E-01 ...
  0.35622E-01  0.11040E+00  0.11176E+00  0.34581E-01  0.11732E+00  0.10628E+00 ...
  0.73667E-01  0.52945E-02  0.63075E-01  0.98566E-01  0.55276E-01  0.19050E-01 ...
  0.10931E+00  0.53875E-01  0.35116E-01];
charge(3,:) = [0.42984E-01  0.97085E-01  0.12359E+00  0.95609E-01  0.46152E-01  0.13519E+00 ...
  0.31062E-01  0.13451E+00  0.13726E+00  0.80789E-01  0.15550E-02  0.78142E-01 ...
  0.10794E+00  0.95881E-01  0.13004E+00  0.16150E-01  0.26977E-02  0.26470E-01 ...
  0.97041E-01  0.17299E-01  0.59700E-01  0.10802E+00  0.32606E-01  0.13035E+00 ...
  0.73704E-01  0.53996E-01  0.10977E+00] ...
  +1i *[0.68190E-01  0.13109E+00  0.13505E+00  0.78894E-01  0.42090E-01  0.25337E-01 ...
  0.62741E-01  0.24619E-01  0.15279E-01  0.32597E-01  0.45932E-01  0.10997E-03 ...
  0.11521E+00  0.99189E-01  0.11171E+00  0.28014E-01  0.24114E-01  0.80821E-01 ...
  0.17738E-01  0.19475E-01  0.12337E+00  0.95289E-01  0.77380E-01  0.11926E+00 ...
  0.12919E+00  0.98115E-01  0.30220E-01];

shift = h/1000;
centers = zeros(3,ns);
for i=1:ns
  centers(1,i) = source(1,i) + shift;
  centers(2,i) = source(2,i);
  centers(3,i) = source(3,i);
end

% now form a multipole expansion at each center
ntm = 7; % what's ntm? Hai
ntot = 0; nterms = zeros(1,nc);
for i = 1:nc
  nterms(i) = ntm + fix(2*cos(pi/2*i));
  ntot = ntot + (nterms(i)+1)*(2*nterms(i)+1); % num of coeffs? 
end

mpole = zeros(1,nd*ntot);
impole(1) = 1;
for i = 1:nc-1
  ilen = (nterms(i)+1)*(2*nterms(i)+1);
  impole(i+1) = impole(i) + nd*ilen;
end 

nlege = 400;
lw = 5*(nlege+1)^2;
%
wlege = zeros(lw,1);
lused = 0;
[wlege, lused] = ylgndrfwini(nlege, wlege, lw, lused);

% keyboard

mpole = zeros(nd*ntot,1); % zinitialize

ns1 = 1;
rscale = 1;
sc = abs(zk)*shift;
if (sc < 1), rscale = sc; end

rscales = zeros(1,nc);
mpole_count = zeros(1,nc);
for i=1:nc-1
  rscales(i) = rscale;
  mpole(impole(i):impole(i+1)-1) = h3dformmpc(nd, zk, rscale, source(:,i), charge(:,i), ns1, centers(:,i), nterms(i), wlege, nlege);  
end
mpole(impole(nc):end) = h3dformmpc(nd, zk, rscale, source(:,nc), charge(:,nc), ns1, centers(:,nc), nterms(nc), wlege, nlege);  

thresh = 1.0d-15;
ifcharge = 1;
ifdipole = 0;
ifpgh = 1;
ntarg = 0;
ifpghtarg = 0;

hess = zeros(nd,6,ns);
targ = zeros(3,ntarg);

keyboard


[pot, grad, hess, pottarg, gradtarg, hesstarg, ier] = hfmm3d(nd, eps, zk, ns, source, ifcharge, ...
      charge, ifdipole, dipvec, [], ifpgh, pot, grad, hess, ntarg, ...
      targ, ifpghtarg);

keyboard
