function [carray,rdpi2,rdmpi2,rdsq3,rdmsq3,dc] = getpwrotmat_mex(nterms,carray,rdpi2,rdmpi2,rdsq3,rdmsq3,dc)
nterms4p1 = 4*nterms+1;
nterms1 = nterms+1;
nterms1s = nterms1^2;
nterms2p1 = nterms*2+1;
rdpi2 = reshape(rdpi2,nterms1s,nterms2p1);
rdmpi2 = reshape(rdmpi2,nterms1s,nterms2p1);
rdsq3 = reshape(rdsq3,nterms1s,nterms2p1);
rdmsq3 = reshape(rdmsq3,nterms1s,nterms2p1);
mex_id_ = 'getpwrotmat(i int64_t[x], io double[xx], io double[xx], io double[xx], io double[xx], io double[xx], io double[xx])';
[carray, rdpi2, rdmpi2, rdsq3, rdmsq3, dc] = laprouts3d(mex_id_, nterms, carray, rdpi2, rdmpi2, rdsq3, rdmsq3, dc, 1, nterms4p1, nterms4p1, nterms1s, nterms2p1, nterms1s, nterms2p1, nterms1s, nterms2p1, nterms1s, nterms2p1, nterms4p1, nterms4p1);
rdpi2 = reshape(rdpi2,nterms1,nterms1,nterms2p1);
rdmpi2 = reshape(rdmpi2,nterms1,nterms1,nterms2p1);
rdsq3 = reshape(rdsq3,nterms1,nterms1,nterms2p1);
rdmsq3 = reshape(rdmsq3,nterms1,nterms1,nterms2p1);
end

