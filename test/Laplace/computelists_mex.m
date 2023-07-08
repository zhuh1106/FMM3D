function [nlist1,list1,nlist2,list2,nlist3,list3,nlist4,list4]=computelists_mex(nlevels,nboxes,laddr,boxsize,centers,iparent,nchild,ichild,isep,nnbors,mnbors,nbors,iper,nlist1,mnlist1,list1,nlist2,mnlist2,list2,nlist3,mnlist3,list3,nlist4,mnlist4,list4)
nlevels1 = nlevels+1;
mex_id_ = 'computelists(i int64_t[x], i int64_t[x], i int64_t[xx], i double[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[xx], i int64_t[x], io int64_t[x], i int64_t[x], io int64_t[xx], io int64_t[x], i int64_t[x], io int64_t[xx], io int64_t[x], i int64_t[x], io int64_t[xx], io int64_t[x], i int64_t[x], io int64_t[xx])';
[nlist1, list1, nlist2, list2, nlist3, list3, nlist4, list4] = laprouts3d(mex_id_, nlevels, nboxes, laddr, boxsize, centers, iparent, nchild, ichild, isep, nnbors, mnbors, nbors, iper, nlist1, mnlist1, list1, nlist2, mnlist2, list2, nlist3, mnlist3, list3, nlist4, mnlist4, list4, 1, 1, 2, nlevels1, nlevels1, 3, nboxes, nboxes, nboxes, 8, nboxes, 1, nboxes, 1, mnbors, nboxes, 1, nboxes, 1, mnlist1, nboxes, nboxes, 1, mnlist2, nboxes, nboxes, 1, mnlist3, nboxes, nboxes, 1, mnlist4, nboxes);
end

