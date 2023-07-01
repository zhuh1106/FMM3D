function  [nlevels,nboxes,ltree] = pts_tree_mem_mex(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree)
if nt == 0 
    targ = zeros(3,1); 
    mex_id_ = 'ptstreemem(i double[xx], i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x])';
[nlevels, nboxes, ltree] = laprouts3d(mex_id_, src, ns, targ, nt, idivflag, ndiv, nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, 3, ns, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
else
    mex_id_ = 'ptstreemem(i double[xx], i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x], io int[x])';
[nlevels, nboxes, ltree] = laprouts3d(mex_id_, src, ns, targ, nt, idivflag, ndiv, nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, 3, ns, 1, 3, nt, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
end % is this ok?

end

