function [cmpole,rmpole,mterms,mpole,impole,ntot,nlege,wlege,lused] = myl3dformmpc(nd, rscale, source, charge, npts, ntm)
% loop over nmpole, call laplace formmpc
%

nmpole = numel(charge)/npts;
cmpole = zeros(3,nmpole);
rmpole = zeros(nmpole,1); % rescaling factor for each mps center
mterms = zeros(nmpole,1); % multipole expansion order 
impole = zeros(nmpole,1); % index of mpole for each mps center
ntot = 0; % total num of multipole expansion coefficients
nlege = 300;
lw = 2*(nlege+1)^2; % dimension
wlege = zeros(lw,1);
lused = 0;
[wlege,lused] = ylgndrfwini(nlege, wlege, lw, lused); % recursion coefficients
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
mpole = zeros(nd*ntot,1); % nd = 1 for now, complex
for i = 1:nmpole
  mpoletmp = zeros(nd,mterms(i)+1,2*mterms(i)+1);
  mpoletmp = l3dformmpc(nd, rscale, source(:,((i-1)*npts+1):i*npts), charge(:,((i-1)*npts+1):i*npts), npts, cmpole(:,i), mterms(i), mpoletmp, wlege, nlege);
%   mpoletmp = reshape(mpoletmp,mterms(i)+1,2*mterms(i)+1);
%   mpoletmp(:,1:mterms(i)) = conj(mpoletmp(:,end:-1:mterms(i)+2)); % use second half multipole expansion
  mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) = mpole(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1) + mpoletmp(:);
end % why imaginary part so small?

end

