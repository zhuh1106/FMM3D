function [U,varargout] = pts_tree3d(srcinfo,varargin)
% comments are from 2d... 07/12/23
%
%   generate level restricted tree based on resolving points
%    (Currently only supports, sorting on sources, sorting on
%      targets, or sorting on sources and targets)
%
%   There is additional functionality to further subdivide
%    an existing tree based on resolving a function, 
%      this functionality tree will be added in later
%
%   This code has the following user callable routines
%
%      pts_tree_mem -> returns memory requirements for creating
%         a tree based on max number of sources/targets
%         in a box (tree length
%         number of boxes, number of levels)
%      pts_tree -> Make the actual tree, returns centers of boxes,
%        colleague info, pts sorted on leaf boxes
%
%   
%  Args:
%
%  -  srcinfo: structure
%        structure containing sourceinfo
%     
%     *  srcinfo.sources: double(2,n)    
%           source locations, $x_{j}$
%
%  Optional args
%  -  targ: double(2,nt)
%        target locations, $t_{i}$ 
%  -  opts: options structure, values in brackets indicate default
%           values wherever applicable
%        opts.ndiv: set number of points for subdivision criterion
%        opts.idivflag: set subdivision criterion (0)
%           opts.idivflag = 0, subdivide on sources only
%           opts.idivflag = 1, subdivide on targets only
%           opts.idivflag = 2, subdivide on sources and targets
%  - ndiv: integer
%        subdivide if relevant number of particles
%        per box is greater than ndiv
%
%  Returns:
%  
%  -  U.nlevels: number of levels
%  -  U.nboxes: number of boxes
%  -  U.ltree: length of tree
%  -  U.itree: integer(ltree) 
%           tree info
%  -  U.iptr: integer(8)
%
%     *  iptr(1) - laddr
%     *  iptr(2) - ilevel
%     *  iptr(3) - iparent
%     *  iptr(4) - nchild
%     *  iptr(5) - ichild
%     *  iptr(6) - ncoll
%     *  iptr(7) - coll
%     *  iptr(8) - ltree
%  -  U.centers: xy coordinates of box centers in the oct tree
%  -  U.boxsize: size of box at each of the levels
%

  src = srcinfo.sources;
  [m,ns] = size(src);
  assert(m==3,'The first dimension of sources must be 3');
  
  if( nargin < 1)
    disp('Not enough input arguments, exiting\n');
    return;
  end
  if( nargin == 1 )
    nt = 0;
    targ = zeros(3,1);
    opts = [];
  elseif (nargin == 2)
    nt = 0;
    targ = zeros(3,1);
    opts = varargin{1};
  elseif (nargin == 3)
    targ = varargin{1};
    [m,nt] = size(targ);
    assert(m==3,'First dimension of targets must be 3');
    opts = [];
  elseif (nargin == 4)
    targ = varargin{1};
    [m,nt] = size(targ);
    assert(m==3,'First dimension of targets must be 3');
    opts = varargin{3};
  end
  ntuse = max(nt,1);

  if ns == 0, disp('Nothing to compute'); return; end;

    
  ndiv = 20;
  idivflag = 0;
  if(isfield(opts,'ndiv'))
    ndiv = opts.ndiv;
  end

  if(isfield(opts,'idivflag'))
    idivflag = opts.idivflag;
  end
  
  nlmin = 0;
  nlmax = 51;
  ifunif = 0; % check results for the case ifunif = 1
  iper = 0;
  nlevels = 0;
  nboxes = 0;
  ltree = 0;

mex_id_ = 'ptstreemem(i double[xx], i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io int64_t[x], io int64_t[x], io int64_t[x])';
[nlevels, nboxes, ltree] = laprouts3d(mex_id_, src, ns, targ, nt, idivflag, ndiv, nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, 3, ns, 1, 3, ntuse, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);

  itree = zeros(1,ltree);
  iptr = zeros(1,8);
  centers = zeros(3,nboxes);
  nlevels1 = nlevels+1;
  boxsize = zeros(1,nlevels1);

mex_id_ = 'ptstreebuild(i double[xx], i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], io int64_t[x], io int64_t[x], io double[xx], io double[x])';
[itree, iptr, centers, boxsize] = laprouts3d(mex_id_, src, ns, targ, nt, idivflag, ndiv, nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, itree, iptr, centers, boxsize, 3, ns, 1, 3, ntuse, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ltree, 8, 3, nboxes, nlevels1);

  n = ns;
  xys = src;
  ixy = zeros(1,n);
  ixyse = zeros(2,nboxes);

mex_id_ = 'ptstreesort(i int64_t[x], i double[xx], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], io int64_t[x], io int64_t[xx])';
[ixy, ixyse] = laprouts3d(mex_id_, n, xys, itree, ltree, nboxes, nlevels, iptr, centers, ixy, ixyse, 1, 3, n, ltree, 1, 1, 1, 8, 3, nboxes, n, 2, nboxes);

  U.nlevels = nlevels;
  U.nboxes = nboxes;
  U.ltree = ltree;
  U.itree = itree;
  U.iptr = iptr;
  U.centers = centers;
  U.boxsize = boxsize;

  varargout{1} = ixy;
  varargout{2} = ixyse;
end