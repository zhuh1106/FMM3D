function arrsort = dreorderf_mex(ndim,n,arr,arrsort,iarr)
len = numel(arr(1,:));
mex_id_ = 'dreorderf(i int64_t[x], i int64_t[x], i double[xx], io double[xx], i int64_t[x])';
[arrsort] = laprouts3d(mex_id_, ndim, n, arr, arrsort, iarr, 1, 1, ndim, len, ndim, len, len);
end

