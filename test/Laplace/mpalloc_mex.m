function [iaddr,lmptot] = mpalloc_mex(nd,laddr,iaddr,nlevels,lmptot,nterms,nboxes)
% one additional imput?
nlevels1 = nlevels+1;
mex_id_ = 'mpalloc0(i int[x], i int[xx], io int[xx], i int[x], io int[x], i int[x], i int[x])';
[iaddr, lmptot] = laprouts3d(mex_id_, nd, laddr, iaddr, nlevels, lmptot, nterms, nboxes, 1, 2, nlevels1, 2, nboxes, 1, 1, nlevels1, 1);
end

