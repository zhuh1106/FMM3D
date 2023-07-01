function [ixy,ixyse] = pts_tree_sort_mex(n,xys,itree,ltree,nboxes,nlevels,iptr,centers,ixy,ixyse)
mex_id_ = 'ptstreesort(i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], io int[x], io int[xx])';
[ixy, ixyse] = laprouts3d(mex_id_, n, xys, itree, ltree, nboxes, nlevels, iptr, centers, ixy, ixyse, 1, 3, n, ltree, 1, 1, 1, 8, 3, nboxes, n, 2, nboxes);
end