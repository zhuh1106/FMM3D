function rmlexp = getpwlistallprocessudnsewexp0_mex(ibox,bs,nboxes,nnbors,nbors,nchild,ichild,centers,isep,nd,ilev,rscale,nterms,iaddrtmp,rmlexp,lmptottmp,rlams,whts,nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,mexpupdown,mexpupphys,mexpdownphys,mexppall,rdplus,rdminus,xs,ys,zs,fexpback,nn,rlsc,rscpow, pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)
ichild=reshape(ichild,8,nboxes);
ndnexptotp=nd*nexptotp;
nboxes6=nboxes*6;
mexp=reshape(mexp,ndnexptotp,nboxes6);
ndnexptotp=nd*nexptotp;
mexppall=reshape(mexppall,ndnexptotp,16);
nterms1 = nterms+1;
nterms1s = nterms1^2;
nterms2p1 = nterms*2+1; 
rdplus = reshape(rdplus,nterms1s,nterms2p1);
rdminus = reshape(rdminus,nterms1s,nterms2p1);
cntlist46 = cntlist4*6;
nw8 = 0;
rlsc = reshape(rlsc,nterms1s,nlams); 
ndnexptot = nd*nexptot;
mexpupdown = reshape(mexpupdown,ndnexptot,2);
% why getpwlistallprocessudnsewexp0 gets Identifier should be a string error is a myth...
%# FORTRAN getpwlistallprocessudnsewexp0(int[1] ibox, double[1] bs, int[1] nboxes, int[1] nnbors, int[nnbors] nbors, int[1] nchild, int[8,nboxes] ichild, double[3,nboxes] centers, int[1] isep, int[1] nd, int[1] ilev, double[1] rscale, int[1] nterms, int[2,nboxes] iaddrtmp, double[lmptottmp] rmlexp, int[1] lmptottmp, double[nlams] rlams, double[nlams] whts, int[1] nlams, int[nlams] nfourier, int[nlams] nphysical, int[1] nthmax, int[1] nexptot, int[1] nexptotp, dcomplex[ndnexptotp,nboxes6] mexp, dcomplex[nd,nexptot] mexpup, dcomplex[nd,nexptot] mexpdown, dcomplex[nd,nexptotp] mexpupphys, dcomplex[nd,nexptotp] mexpdownphys, dcomplex[ndnexptotp,16] mexppall, double[nterms1s,nterms2p1] rdplus, double[nterms1s,nterms2p1] rdminus, dcomplex[11,nexptotp] xs, dcomplex[11,nexptotp] ys, double[5,nexptotp] zs, dcomplex[nn] fexpback, int[1] nn, double[nterms1s,nlams] rlsc, double[nterms1] rscpow, dcomplex[ndnexptotp,cntlist46] pgboxwexp, int[1] cntlist4, int[nboxes] list4ct, int[nboxes] nlist4, int[mnlist4,nboxes] list4, int[1] mnlist4);
mex_id_ = 'getpwlistallprocessudnsewexp02(i int64_t[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[xx], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], i int64_t[xx], io double[x], i int64_t[x], i double[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i double[x], i dcomplex[x], i int64_t[x])';
[rmlexp] = laprouts3d(mex_id_, ibox, bs, nboxes, nnbors, nbors, nchild, ichild, centers, isep, nd, ilev, rscale, nterms, iaddrtmp, rmlexp, lmptottmp, rlams, whts, nlams, nfourier, nphysical, nthmax, nexptot, nexptotp, mexp, mexpupdown, mexpupphys, mexpdownphys, mexppall, rdplus, rdminus, xs, ys, zs, rlsc, rscpow, fexpback, nn, 1, 1, 1, 1, nnbors, 1, 8, nboxes, 3, nboxes, 1, 1, 1, 1, 1, 2, nboxes, lmptottmp, 1, nlams, nlams, 1, nlams, nlams, 1, 1, 1, ndnexptotp, nboxes6, ndnexptot, 2, nd, nexptotp, nd, nexptotp, ndnexptotp, 16, nterms1s, nterms2p1, nterms1s, nterms2p1, 11, nexptotp, 11, nexptotp, 5, nexptotp, nterms1s, nlams, nterms1, nn, 1);
end
