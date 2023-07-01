function [pot,grad,hess,pottarg,gradtarg,hesstarg] = lfmm3d_mex(nd,eps,ns,source,ifcharge,charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,nt,targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
% need to take care of ns or nt = 0 ?
nd3 = nd*3;
nd6 = nd*6;
dipvec = reshape(dipvec,nd3,ns);
grad = reshape(grad,nd3,ns);
hess = reshape(hess,nd6,ns);
gradtarg = reshape(gradtarg,nd3,nt);
hesstarg = reshape(hesstarg,nd6,nt);
mex_id_ = 'lfmm3d(i int[x], i double[x], i int[x], i double[xx], i int[x], i double[xx], i int[x], i double[xx], i int[x], i int[x], io double[xx], io double[xx], io double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i int[x])';
[pot, grad, hess, pottarg, gradtarg, hesstarg] = laprouts3d(mex_id_, nd, eps, ns, source, ifcharge, charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, nt, targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier, 1, 1, 1, 3, ns, 1, nd, ns, 1, nd3, ns, 1, 1, nd, ns, nd3, ns, nd6, ns, 1, 3, nt, 1, nd, nt, nd3, nt, nd6, nt, 1);
grad = reshape(grad,nd,3,ns);
hess = reshape(hess,nd,6,ns);
gradtarg = reshape(gradtarg,nd,3,nt);
hesstarg = reshape(hesstarg,nd,6,nt);
end

