function scarray = l3dmpevalhessdini_mex(nterms,scarray)
nscarray = length(scarray);
mex_id_ = 'l3dmpevalhessdini(i int64_t[x], io double[x])';
[scarray] = laprouts3d(mex_id_, nterms, scarray, 1, nscarray);
end

