function dc = getsqrtbinomialcoeffs_mex(n,dc)
n1 = n+1;
mex_id_ = 'getsqrtbinomialcoeffs(i int[x], io double[xx])';
[dc] = laprouts3d(mex_id_, n, dc, 1, n1, n1);
end

