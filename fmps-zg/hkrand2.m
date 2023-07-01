function x = hkrand2(iseed_hk)
%
%
mex_id_ = 'hkrand2(i int[x], o double[x])';
[x] = fmps(mex_id_, iseed_hk, 1, 1);

