function [itree,iptr,centers,boxsize] = pts_tree_build_mex(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree,itree,iptr,centers,boxsize)
nlevels1 = nlevels+1;
ntuse = max(nt,1);
if nt == 0, targ = zeros(3,1); end
mex_id_ = 'ptstreebuild(i double[xx], i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io int64_t[x], io int64_t[x], io double[xx], io double[x])';
[itree, iptr, centers, boxsize] = laprouts3d(mex_id_, src, ns, targ, nt, idivflag, ndiv, nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, itree, iptr, centers, boxsize, 3, ns, 1, 3, ntuse, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ltree, 8, 3, nboxes, nlevels1);
%ptstreebuild(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree,itree,iptr,centers,boxsize);
end

