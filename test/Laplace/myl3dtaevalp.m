function pot = myl3dtaevalp(nd, rmpole, cmpole, impole, local, mterms, source, npts, wlege, nlege)
%
ns = numel(source(1,:));
nmpole = ns/npts;
pot = zeros(nd,ns);
for i = 1:nmpole % loop over mps center, ship local expansion to targ
  poti = zeros(nd,npts);
  locali = reshape(local(impole(i):impole(i)+nd*(mterms(i)+1)*(2*mterms(i)+1)-1),mterms(i)+1,2*mterms(i)+1);
%   locali(:,1:mterms(i)) = 0; 
%   locali(locali==0) = 1; % ... gmres initialze things to 0, this is ot a godd way to assign 1/2 of coeffs to 0
  poti = l3dtaevalp(nd, rmpole(i), cmpole(:,i), locali, mterms(i), source(:,((i-1)*npts+1):i*npts), npts, poti, wlege, nlege);
  pot(:,((i-1)*npts+1):i*npts) = pot(:,((i-1)*npts+1):i*npts) + poti;
end
end