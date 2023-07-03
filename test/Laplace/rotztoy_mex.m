function mrotate = rotztoy_mex(nd,nterms,mpole,mrotate,rdminus)
ndnterms = nd*(nterms+1);
ntermst2p1 = 2*nterms+1;
mpole = reshape(mpole,ndnterms,ntermst2p1);
mrotate = reshape(mrotate,ndnterms,ntermst2p1);
nterms1 = nterms+1;
nterms1s = nterms1^2;
nterms2p1 = nterms*2+1; 
rdminus = reshape(rdminus,nterms1s,nterms2p1);
mex_id_ = 'rotztoy(i int[x], i int[x], i dcomplex[xx], io dcomplex[xx], i double[xx])';
[mrotate] = laprouts3d(mex_id_, nd, nterms, mpole, mrotate, rdminus, 1, 1, ndnterms, ntermst2p1, ndnterms, ntermst2p1, nterms1s, nterms2p1);
mrotate = reshape(mrotate,nd,nterms+1,2*nterms+1);
end

