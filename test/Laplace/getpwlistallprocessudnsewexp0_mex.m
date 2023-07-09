function rmlexp = getpwlistallprocessudnsewexp0_mex(ibox,bs,nboxes,nnbors,nbors,nchild,ichild,centers,isep,nd,ilev,rscale,nterms,iaddrtmp,rmlexp,lmptottmp,rlams,whts,nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,mexpup,mexpdown,mexpupphys,mexpdownphys,mexppall,rdplus,rdminus,xs,ys,zs,fexpback,nn,rlsc,rscpow, pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)
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
pgboxwexp = reshape(pgboxwexp,ndnexptotp,cntlist46);
mex_id_ = 'getpwlistallprocessudnsewexp0(i int64_t[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[xx], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], i int64_t[xx], io double[x], i int64_t[x], i double[x], i double[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i dcomplex[xx], i double[xx], i double[xx], i dcomplex[xx], i dcomplex[xx], i double[xx], i dcomplex[x], i int64_t[x], i double[xx], i double[x], i dcomplex[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[xx], i int64_t[x])';
[rmlexp] = laprouts3d(mex_id_, ibox, bs, nboxes, nnbors, nbors, nchild, ichild, centers, isep, nd, ilev, rscale, nterms, iaddrtmp, rmlexp, lmptottmp, rlams, whts, nlams, nfourier, nphysical, nthmax, nexptot, nexptotp, mexp, mexpup, mexpdown, mexpupphys, mexpdownphys, mexppall, rdplus, rdminus, xs, ys, zs, fexpback, nn, rlsc, rscpow, pgboxwexp, cntlist4, list4ct, nlist4, list4, mnlist4, 1, 1, 1, 1, nnbors, 1, 8, nboxes, 3, nboxes, 1, 1, 1, 1, 1, 2, nboxes, lmptottmp, 1, nlams, nlams, 1, nlams, nlams, 1, 1, 1, ndnexptotp, nboxes6, nd, nexptot, nd, nexptot, nd, nexptotp, nd, nexptotp, ndnexptotp, 16, nterms1s, nterms2p1, nterms1s, nterms2p1, 11, nexptotp, 11, nexptotp, 5, nexptotp, nn, 1, nterms1s, nlams, nterms1, ndnexptotp, cntlist46, 1, nboxes, nboxes, mnlist4, nboxes, 1);
end

