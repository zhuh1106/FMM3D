function rlsc = rlscini_mex(rlsc,nlambs,rlams,nterms)
nterms1 = nterms+1;
nterms1s = nterms1^2;
rlsc = reshape(rlsc,nterms1s,nlambs);
mex_id_ = 'rlscini(io double[xx], i int[x], i double[x], i int[x])';
[rlsc] = laprouts3d(mex_id_, rlsc, nlambs, rlams, nterms, nterms1s, nlambs, 1, nlambs, 1);
rlsc = reshape(rlsc,nterms1,nterms1,nlambs);
end

