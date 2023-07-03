function arrsort = ireorderf_mex(ndim,n,arr,arrsort,iarr)
len = numel(arr(1,:));
mex_id_ = 'ireorderf(i int[x], i int[x], i int[xx], io int[xx], i int[x])';
[arrsort] = laprouts3d(mex_id_, ndim, n, arr, arrsort, iarr, 1, 1, ndim, len, ndim, len, len);
end

