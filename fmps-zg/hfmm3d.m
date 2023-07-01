function [pot, grad, hess, pottarg, gradtarg, hesstarg, ier] = hfmm3d(nd,eps,zk,ns,source,ifcharge,charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg)
%
%
[m,nt] = size(targ);
iper = 0;
nd3 = 3*nd;
nd6 = 6*nd;
ntuse = max(nt,1);
if nt == 0, targ = zeros(3,1); end
grad = reshape(grad,[nd3,ns]);
hess = reshape(hess,[nd6,ns]);
dipvec = reshape(dipvec,[nd3,ns]);
pottarg = complex(zeros(nd,ntuse));
gradtarg = complex(zeros(nd3,ntuse));
hesstarg = complex(zeros(nd6,ntuse));
ier = 0;
mex_id_ = 'hfmm3d(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], io dcomplex[xx], io int[x])';
[pot, grad, hess, pottarg, gradtarg, hesstarg, ier] = fmps(mex_id_, nd, eps, zk, ns, source, ifcharge, charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, ntuse, targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier, 1, 1, 1, 1, 3, ns, 1, nd, ns, 1, nd3, ns, 1, 1, nd, ns, nd3, ns, nd6, ns, 1, 3, ntuse, 1, nd, ntuse, nd3, ntuse, nd6, ntuse, 1);
% hfmm3d(int[1] nd, double[1] eps, dcomplex[1] zk, int[1] ns, double[3,ns] source, int[1] ifcharge, dcomplex[nd,ns] charge, int[1] ifdipole, dcomplex[nd3,ns] dipvec, int[1] iper, int[1] ifpgh, inout dcomplex[nd,ns] pot, inout dcomplex[nd3,ns] grad, inout dcomplex[nd6,ns] hess, int[1] ntarg, targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)

