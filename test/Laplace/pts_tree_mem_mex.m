function  [nlevels,nboxes,ltree] = pts_tree_mem_mex(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree)
ntuse = max(nt,1);
if nt == 0, targ = zeros(3,1); end
mex_id_ = 'ptstreemem(i double[xx], i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io int64_t[x], io int64_t[x], io int64_t[x])';
[nlevels, nboxes, ltree] = laprouts3d(mex_id_, src, ns, targ, nt, idivflag, ndiv, nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, 3, ns, 1, 3, ntuse, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
end

