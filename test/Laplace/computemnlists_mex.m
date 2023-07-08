function [mnlist1,mnlist2,mnlist3,mnlist4] = computemnlists_mex(nlevels,nboxes,laddr,boxsize,centers,iparent,nchild,ichild,isep,nnbors,mnbors,nbors,iper,mnlist1,mnlist2,mnlist3,mnlist4)
nlevels1 = nlevels+1;
mex_id_ = 'computemnlists(i int64_t[x], i int64_t[x], i int64_t[xx], i double[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[xx], i int64_t[x], io int64_t[x], io int64_t[x], io int64_t[x], io int64_t[x])';
[mnlist1, mnlist2, mnlist3, mnlist4] = laprouts3d(mex_id_, nlevels, nboxes, laddr, boxsize, centers, iparent, nchild, ichild, isep, nnbors, mnbors, nbors, iper, mnlist1, mnlist2, mnlist3, mnlist4, 1, 1, 2, nlevels1, nlevels1, 3, nboxes, nboxes, nboxes, 8, nboxes, 1, nboxes, 1, mnbors, nboxes, 1, 1, 1, 1, 1);
end

