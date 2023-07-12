function rmlexp = getlist3pwallprocessudnsewexp0_mex(ibox,bs,nboxes,nlist3,list3,isep,centers,nd,nterms,rmlexp,rlams,whts,nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,mexpup,mexpdown,mexpupphys,mexpdownphys,mexppall1,mexppall2,xs,ys,zs,fexpback,nn,rlsc,rscpow,rdplus,rdminus)
nboxes6=nboxes*6;
ndnexptotp=nd*nexptotp;
mexp=reshape(mexp,ndnexptotp,nboxes6);
ndnterms1t2nterms1 = nd*(nterms+1)*(2*nterms+1);
nterms1 = nterms+1;
nterms1s = nterms1^2;
nterms2p1 = nterms*2+1; 
rdplus = reshape(rdplus,nterms1s,nterms2p1);
rdminus = reshape(rdminus,nterms1s,nterms2p1);
rlsc = reshape(rlsc,nterms1s,nlams); 
mex_id_ = 'getlist3pwallprocessudnsewexp0(i int64_t[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], io dcomplex[xx], i double[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i double[xx], i dcomplex[x], i int64_t[x], i double[xx], i double[x], i double[xx], i double[xx])';
[rmlexp] = laprouts3d(mex_id_, ibox, bs, nboxes, nlist3, list3, isep, centers, nd, nterms, rmlexp, rlams, whts, nlams, nfourier, nphysical, nthmax, nexptot, nexptotp, mexp, mexpup, mexpdown, mexpupphys, mexpdownphys, mexppall1, mexppall2, xs, ys, zs, fexpback, nn, rlsc, rscpow, rdplus, rdminus, 1, 1, 1, 1, nlist3, 1, 3, nboxes, 1, 1, ndnterms1t2nterms1, 8, nlams, nlams, 1, nlams, nlams, 1, 1, 1, ndnexptotp, nboxes6, nd, nexptot, nd, nexptot, nd, nexptotp, nd, nexptotp, nd, nexptotp, nd, nexptotp, 11, nexptotp, 11, nexptotp, 5, nexptotp, nn, 1, nterms1s, nlams, nterms1, nterms1s, nterms2p1, nterms1s, nterms2p1);
% getlist3pwallprocessudnsewexp0(               ibox,           bs,        nboxes,        nlist3,             list3,        isep,                  centers,        nd,        nterms,                                      rmlexp,               rlams,               whts,        nlams,            nfourier,            nphysical,        nthmax,        nexptot,        nexptotp,                              mexp,                      mexpup,                      mexpdown,                       mexpupphys,                       mexpdownphys,                       mexppall1,                       mexppall2,                       xs,                       ys,                    zs,              fexpback,        nn,                        rlsc,                 rscpow,                            rdplus,                            rdminus);
rmlexp; % reshape?
end

