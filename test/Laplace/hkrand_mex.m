function hkrand = hkrand_mex(iseed_hk)
mex_id_ = 'o double = hkrand(i int[x])';
[hkrand] = laprouts3d(mex_id_, iseed_hk, 1);
end

