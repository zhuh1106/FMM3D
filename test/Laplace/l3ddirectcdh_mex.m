function [pottarg, gradtarg, hesstarg] = l3ddirectcdh_mex(nd, sources, charges, dipoles, ns, targ, nt, pottarg, gradtarg, hesstarg, thresh)
nd3 = 3*nd;
nd6 = 6*nd;
dipoles = reshape(dipoles,nd3,ns);
gradtarg = reshape(gradtarg,nd3,nt);
hesstarg = reshape(hesstarg,nd6,nt);
mex_id_ = 'l3ddirectcdh(i int[x], i double[xx], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = laprouts3d(mex_id_, nd, sources, charges, dipoles, ns, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, nd6, nt, 1);
gradtarg = reshape(gradtarg,nd,3,nt);
hesstarg = reshape(hesstarg,nd,6,nt);
end

