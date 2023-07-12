clear srcinfo

ns = 4000;

idist=1;
srcinfo.sources = rand(3,ns);

if( idist == 1 )
theta=rand(1,ns)*pi;
phi=rand(1,ns)*2*pi;
srcinfo.sources(1,:)=.5*cos(phi).*sin(theta);
srcinfo.sources(2,:)=.5*sin(phi).*sin(theta);
srcinfo.sources(3,:)=.5*cos(theta);
end

srcinfo.charges = rand(1,ns);

nt = 3999;
targ = rand(3,nt);

eps = 1e-5;
ntests = 54;
ipass = zeros(ntests,1);
errs = zeros(ntests,1);
ntest = 10;

stmp = srcinfo.sources(:,1:ntest);
ttmp = targ(:,1:ntest);

itest = 1;

% Test tree
opts.ndiv = 40;
[U,ixy,ixyse] = pts_tree3d(srcinfo,opts);

  nlevels = U.nlevels;
  nboxes = U.nboxes;
  ltree = U.ltree;
  itree = U.itree;
  iptr = U.iptr;
  centers = U.centers;
  boxsize = U.boxsize;

% plot
  level = itree(iptr(2):iptr(3)-1);
  nchild = itree(iptr(4):iptr(5)-1);
  sources_sort = srcinfo.sources(:,ixy);
  figure(1),clf,
  for k=1:nboxes
    if nchild(k)==0 % no children
      nodesXud = centers(1,k) + boxsize(level(k))/2*[[1,-1,-1,1,1] [1,-1,-1,1,1]]/2; % why divided by 4?
      nodesYud = centers(2,k) + boxsize(level(k))/2*[[1,1,-1,-1,1] [1,1,-1,-1,1]]/2;
      nodesZud = centers(3,k) + boxsize(level(k))/2*[[-1,-1,-1,-1,-1] [1,1,1,1,1]]/2;

      nodesXew = centers(1,k) + boxsize(level(k))/2*[[-1,-1,-1,-1,-1] [1,1,1,1,1]]/2; % why divided by 4?
      nodesYew = centers(2,k) + boxsize(level(k))/2*[[1,1,-1,-1,1] [1,1,-1,-1,1]]/2;
      nodesZew = centers(3,k) + boxsize(level(k))/2*[[1,-1,-1,1,1] [1,-1,-1,1,1]]/2;
      %
      if ixyse(2,k) > ixyse(1,k) % plot pts box by box (some does not have any)
        sources_sort_k = sources_sort(:,ixyse(1,k):ixyse(2,k));
        if 0
          plot3(sources_sort_k(1,:),sources_sort_k(2,:),sources_sort_k(3,:),'.'); hold on  
        else
          plot3(nodesXud,nodesYud,nodesZud,'-k'); hold on
          plot3(nodesXew,nodesYew,nodesZew,'-k');
          plot3(centers(1,k),centers(2,k),centers(3,k),'*r'), hold on,
        end
        pause(0.01)
      end
%       keyboard
    end
  end
  plot3(srcinfo.sources(1,:),srcinfo.sources(2,:),srcinfo.sources(3,:),'.b')
  title('plot of leaf boxes')
  axis equal
keyboard