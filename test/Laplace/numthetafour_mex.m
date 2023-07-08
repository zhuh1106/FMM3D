function numtets = numthetafour_mex(numtets,nlams)
mex_id_ = 'numthetafour(io int64_t[x], i int64_t[x])';
[numtets] = laprouts3d(mex_id_, numtets, nlams, nlams, 1);
end

