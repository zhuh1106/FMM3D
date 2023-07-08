function [iaddr,lmptot] = mpalloc_mex(nd,laddr,iaddr,nlevels,lmptot,nterms,nboxes)
% one additional imput?
nlevels1 = nlevels+1;
mex_id_ = 'mpalloc0(i int64_t[x], i int64_t[xx], io int64_t[xx], i int64_t[x], io int64_t[x], i int64_t[x], i int64_t[x])';
[iaddr, lmptot] = laprouts3d(mex_id_, nd, laddr, iaddr, nlevels, lmptot, nterms, nboxes, 1, 2, nlevels1, 2, nboxes, 1, 1, nlevels1, 1);
end

