function [ixy,ixyse] = pts_tree_sort_mex(n,xys,itree,ltree,nboxes,nlevels,iptr,centers,ixy,ixyse)
mex_id_ = 'ptstreesort(i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], io int64_t[x], io int64_t[xx])';
[ixy, ixyse] = laprouts3d(mex_id_, n, xys, itree, ltree, nboxes, nlevels, iptr, centers, ixy, ixyse, 1, 3, n, ltree, 1, 1, 1, 8, 3, nboxes, n, 2, nboxes);
end

