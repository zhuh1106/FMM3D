function [mexpupf,mexpdownf] = mpoletoexp_mex(nd,mpole,nterms,nlambs,numtets,nexptot,mexpupf,mexpdownf,rlsc)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
nterms1 = nterms+1;
nterms1s = nterms1^2;
rlsc = reshape(rlsc,nterms1s,nlambs);
mex_id_ = 'mpoletoexp(i int64_t[x], i dcomplex[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io dcomplex[xx], io dcomplex[xx], i double[xx])';
[mexpupf, mexpdownf] = laprouts3d(mex_id_, nd, mpole, nterms, nlambs, numtets, nexptot, mexpupf, mexpdownf, rlsc, 1, ndnterms, ntermst2p1, 1, 1, nlambs, 1, nd, nexptot, nd, nexptot, nterms1s, nlambs);
end

