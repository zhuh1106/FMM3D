function [u,v] = gl_grid(p)
% GL_GRID - Returns the Gauss-Legendre--uniform grid on the unit sphere
% [0,pi]x[0,2pi).
%
% SEE ALSO: G_GRID
%
  
[nu,nv]  = spharm_grid_size(p);

lambda = (0:nv-1)'*2*pi/nv;
theta  = acos(g_grid(nu));
[v, u] = meshgrid(lambda,theta);
u = u(:); v = v(:);
end

function [nu,nv, p]=spharm_grid_size(p,ntot)
% SPHARM_GRID_SIZE returns the lattitude and longitude grid size for
% spherical harmonic order p.

if(nargin<2), ntot=-1;end
if(isempty(p))
    p  = (round(sqrt(2*ntot+1))-1)/2;
end
nu = p + 1;
nv = 2*p;

if(nargin>1), assert(nv*nu==ntot);end
end