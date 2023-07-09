function [isorted,iboxfl,subcenters] = subdividebox_mex(pos,npts,center,boxsize,isorted,iboxfl,subcenters)
nmaxt = length(isorted);
mex_id_ = 'subdividebox(i double[xx], i int64_t[x], i double[x], i double[x], io int64_t[x], io int64_t[xx], io double[xx])';
[isorted, iboxfl, subcenters] = laprouts3d(mex_id_, pos, npts, center, boxsize, isorted, iboxfl, subcenters, 3, npts, 1, 3, 1, nmaxt, 2, 8, 3, 8);
end

