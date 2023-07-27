function [S,S_up] = mySurfaceSph(p)
S = mySurfaceSph0(p);
S_up = mySurfaceSph0(p+1);
end

function S = mySurfaceSph0(p)
% Shape = 'sphere'
[u, v] = gl_grid(p);
ax = 1; ay = 1; az = 1; rho = 1;
x = ax*rho.*sin(u).*cos(v);
y = ay*rho.*sin(u).*sin(v);
z = az*rho.*cos(u);
S = [];
S.cart.x = x; S.cart.y = y; S.cart.z = z;
S.upFreqCoeff = 1;
S.filterFreqCoeff = 1;
S.stokesOpStale = 1;
S.dblLayerOpStale = 1;
S.p = p;
S.shc.x = shAna(S.cart.x); S.shc.y = shAna(S.cart.y); S.shc.z = shAna(S.cart.z);
S.upFreq = p;
S.filterFreq = p;
[gp,gpAll]= calcGeoProp(S);
S.geoProp = gp;
% a few additional 
S.x = [S.cart.x,S.cart.y,S.cart.z]'; % source on the boundary surface
S.nx = [S.geoProp.nor.x S.geoProp.nor.y S.geoProp.nor.z]';
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
S.w = (S.geoProp.W.*wt(:))'; % quadr weights (GL x Trapezoidal)
end

function [gp,gpAll]= calcGeoProp(S)
% not all variables are checked
%
% gp = calcGeoProp(X)  calculates the geometric properties of the srurface
% defined by X. X is of class vesicle. gp is a structure containing the
% first fundamental coefficients E,F,G, second fundamental coefficients
% L,M,N, surface normal vector nor, area element W, mean and Gaussian
% curvatures H,K, spherical Laplacian Lap, as well as the function handle
% for surface Gradient and surface divergence.
%

  % NOTES: Bellow, gpAll holds the geoProps in the upsampled frequency. The
  % fundamental forms are not continuous and cannot be interpolated via
  % spherical harmonics, it is more accurate to calculate them in the original
  % shOrder. The geoProps in the same spherical order of S (or X) are stored in
  % gp.
  %
  % Also note that interpolating the normal causes it to lose its unity size
  % (n.n=1); if you decide to interpolate, you need to renormalize.
  %

  %-- Function handles for filtering and upsampling depending on the options
  % fH is the filtering handle and 
  % dH is the differentiation handle 
  % X is upsampled based on options
  if(S.upFreq>S.p && S.filterFreq==S.p)
    X = interpsh(S.cart, S.upFreq); 
    dH = @(f,dv,du) DmSH(interpsh(f, S.upFreq), dv, du);
    fH = @(f) interpsh(f, S.p);
  elseif(S.upFreq>S.p)
    X = interpsh(S.cart, S.upFreq); 
    dH = @(f,dv,du) DmSH(interpsh(f, S.upFreq), dv, du);
    fH = @(f) filtersh(interpsh(f,S.p), S.filterFreq);
  else
    X = S.cart;
    dH = @DmSH;
    fH = @(f) filtersh(f, S.filterFreq);
  end
  
  %-- Taking derivatives (in pUp)
  gpAll.Xu  = myDmSH(X,0,1); gpAll.Xv  = myDmSH(X,1,0);
  gpAll.Xuu = myDmSH(X,0,2); gpAll.Xuv = myDmSH(X,1,1); gpAll.Xvv = myDmSH(X,2,0);

  %-- First fundamental coeffs
  gpAll.E = mydot(gpAll.Xu, gpAll.Xu); 
  gpAll.F = mydot(gpAll.Xu, gpAll.Xv); 
  gpAll.G = mydot(gpAll.Xv, gpAll.Xv);
  W2 = abs(gpAll.E.*gpAll.G - gpAll.F.^2); gpAll.W = sqrt(W2);
  gpAll.nor = mycross(gpAll.Xu,gpAll.Xv); gpAll.nor = myrdivide(gpAll.nor,gpAll.W);
  
  %-- Second fundamental coefficients, not generally continuous so cannot be
  %interpolated down
  gpAll.L = mydot(gpAll.Xuu, gpAll.nor); 
  gpAll.M = mydot(gpAll.Xuv, gpAll.nor);
  gpAll.N = mydot(gpAll.Xvv, gpAll.nor);
  
  %-- Gauss and mean curvatures
  gpAll.H = 0.5*(gpAll.E.*gpAll.N - 2*gpAll.F.*gpAll.M + gpAll.G.*gpAll.L)./W2;
  gpAll.K = (gpAll.L.*gpAll.N - gpAll.M.*gpAll.M)./W2;

  %-- Mapping back to p 
  gp.H = fH(gpAll.H); gp.K = fH(gpAll.K);

  %-- Calculating the discontinuous functions in p frequency
  gp.Xu = myDmSH(S.cart,0,1);
  gp.Xv = myDmSH(S.cart,1,0);
  % normal is continuous, but interpolating it will cause it to loose its unit
  % length. Either we can calculate in p or renormalize after downsampling. We
  % calculate it in p, because it is Xu and Xv are spectrally accurate.
  gp.nor = mycross(gp.Xu, gp.Xv); 
  gp.W = mynorm(gp.nor);
  gp.nor = myrdivide(gp.nor,gp.W);
  
  %-- Surface operators coefficients
  Cu = myrdivide(myminus(mytimes(gpAll.G,gpAll.Xu),mytimes(gpAll.F,gpAll.Xv)),W2);
  Cv = myrdivide(myminus(mytimes(gpAll.E,gpAll.Xv),mytimes(gpAll.F,gpAll.Xu)),W2);

  gp.Grad_noisy = @(f) dH(f,0,1).*Cu + dH(f,1,0).*Cv;
  gp.Div_noisy  = @(F) sum(dH(F,0,1).*Cu + dH(F,1,0).*Cv,2);

  %-- Surface operators
  gp.Grad = @(f) fH(gp.Grad_noisy(f));
  gp.Div  = @(F) fH(gp.Div_noisy(F));
  gp.Laplacian  = @(f) fH(gp.Div_noisy(gp.Grad_noisy(f)));

  %-- Linearized curvature
  gp.linCurv = @(X) fH(( gpAll.E.*dot(dH(X,2,0),gpAll.nor) - ...
                       2*gpAll.F.*dot(dH(X,1,1),gpAll.nor) + ...
                         gpAll.G.*dot(dH(X,0,2),gpAll.nor))./W2/2);

  %-- Linearized bending and tension operators on the surface
  H2K = 2*(gpAll.H.^2-gpAll.K);
  gp.bendingCoeff = @(f) -(gp.Laplacian(f) + fH(interpsh(f, S.upFreq).*H2K));
  gp.bendingOp = @(f) gp.bendingCoeff(f).*gp.nor;
  gp.tensionOp = @(f) gp.Grad(f) + 2*fH(interpsh(f,S.upFreq).*gpAll.H).*gp.nor;

  %Ordering fields 
  gp = orderfields(gp);
end

function vec = myDmSH(vec,dv,du)
  c = DmSH([vec.x vec.y vec.z],dv,du);
  vec.x = c(:,1);
  vec.y = c(:,2);
  vec.z = c(:,3);
end

function nv = mynorm(vec)
  nv = sqrt(mydot(vec,vec));
end

function vec1 = myminus(vec1,vec2)
  vec1.x = vec1.x-vec2.x;
  vec1.y = vec1.y-vec2.y;
  vec1.z = vec1.z-vec2.z;
end

function val = mydot(vec1,vec2)
  val = vec1.x.*conj(vec2.x) + vec1.y.*conj(vec2.y) + vec1.z.*conj(vec2.z);
end

function vec = mytimes(scalar,vec)
  %-- elementwise multioication .*
  if(isa(scalar,'vec3d'))
    vec.x = scalar.x.*vec.x;
    vec.y = scalar.y.*vec.y;
    vec.z = scalar.z.*vec.z;
  else
    vec.x = scalar.*vec.x;
    vec.y = scalar.*vec.y;
    vec.z = scalar.*vec.z;
  end
end

function val = mycross(vec1,vec2)
  val.x =  vec1.y.*vec2.z - vec1.z.*vec2.y;
  val.y = -vec1.x.*vec2.z + vec1.z.*vec2.x;
  val.z =  vec1.x.*vec2.y - vec1.y.*vec2.x;
end

function vec = myrdivide(vec,scalar)
  %-- elementwise division ./
  if(isa(scalar,'vec3d'))
    vec.x = vec.x./scalar.x;
    vec.y = vec.y./scalar.y;
    vec.z = vec.z./scalar.z;
  else
    vec.x = vec.x./scalar;
    vec.y = vec.y./scalar;
    vec.z = vec.z./scalar;
  end
end

function dF = DmSH(F,dv,du,useMats)
% DMSH - Takes the derivative of function(s) defined on the sphere with
% respect to parameters through spherical harmonic transformation. Each
% column of F is a function defined on a Gauss-Legendre-uniform grid on the
% sphere. The grid may be produced by PARDOMAIN().
%
% DmSH(F,dv,du) : takes the derivative of the function F in the azimuth
% direction(east-west) dv times and in the elevation(south-north)
% direction du times (du<2).
%
% DmSH(F,dv,du, useMats) : By default DmSH uses matrices for
% differentation, if they are not present, it generates and saves them for
% future use to the '/data/' folder. However, this may slow down the first
% evaluation of the function for a given size. Set useMats to false for
% faster evaluation for large vectors.
%
% SEE ALSO: PARDOMAIN, SHANA, SHSYN.
%

  persistent Du Duu Dv freq

  pMax = 49;
  mem = whos('Du');
  if(mem.bytes > 4e8), Du=[]; Duu=[]; Dv=[]; freq = []; end

  %-- Checking the input
  if(nargin==0), testDmSH(); return;end
  if(du>2)
    error(['Only first and second order differentiation is permitted' ...
           ' for elevation']);
  end

  %-- Getting the size
  [d1, d2] = size(F);
  p = (sqrt(2*d1+1)-1)/2;
  if(p~=fix(p)),
    error(['The input function should be defined on a' ...
           ' Gauss-Legendre--uniform grid. Check the size of input.']);
  end
  if(nargin<4), useMats = true;end

  if(useMats && p<pMax)
    idx = find(freq==p);
    if(isempty(idx))
      [M1, M2, M3] = getMats(p);
      Du{end+1} = M1;
      Dv{end+1} = M2;
      Duu{end+1} = M3;
      freq(end+1) = p;
      idx = length(freq);
    end

    dF = F;
    for k=1:dv, dF = Dv{idx}*dF;end
    if(du==1)
      dF = Du{idx}*dF;
    elseif(du==2)
      dF = Duu{idx}*dF;
    end

  else
    printMsg('* Direct differentiation via SHT for p=%d and (du,dv)=(%d,%d).\n',p,du,dv);
    SHC = shAna(F);
    if(dv>0)
      mMat = (1i*repmat(-p:p,p+1,1)).^dv;
      mMat = shrinkShVec(mMat(:));
      SHC = repmat(mMat,1,d2).*SHC;
    end

    switch du
     case 0
      dF = shSyn(SHC);
     case 1
      dF = sumBasis(SHC, 'Hnm');
     case 2
      dF = sumBasis(SHC, 'Wnm');
    end
    printMsg('* Direct differentiation Successful\n');
  end

  if( isreal( F ) ), dF = real( dF ); end
end

function [Du, Dv, Duu] = getMats(p)

  np = 2*p*(p+1);
  printMsg('* Reading the differentiation matrices for p=%d:',p);
  [Du,  m1] = readData('Du'  ,p,[np np]);
  [Dv,  m2] = readData('Dv'  ,p,[np np]);
  [Duu, m3] = readData('Duu' ,p,[np np]);

  if(m1 && m2 && m3)
    printMsg(' Successful\n');
    return;
  end

  printMsg(' Failed\n');
  printMsg('* Generating the differentiation matrices for p=%d:\n', p);

  useMat = false; I = eye(np);
  if(~m1), Du  = DmSH(I, 0, 1, useMat); writeData('Du' ,p,Du );end
  if(~m2), Dv  = DmSH(I, 1, 0, useMat); writeData('Dv' ,p,Dv );end
  if(~m3), Duu = DmSH(I, 0, 2, useMat); writeData('Duu',p,Duu);end
end

function [ff,FmatRet] = filtersh(f,freq,useMat)
% FILTERSH(f,freq) - filters the function f such that the highest spherical
% harmonic frequency that is present after filtering is freq. Therefore,
% filtersh(f,p) where p is the current band limit of f does not modify f.
%

  persistent Fmat matFreqs

  mem = whos('Fmat');
  if(mem.bytes > 4e8), Fmat = []; matFreqs = []; end
  
  %%-- Checking the input
  if(nargin==0), testFilterSh(); return; end
  np = size(f,1); p = (sqrt(2*np+1)-1)/2;
  
  if(nargin<2 || isempty(freq)), freq = p;end
  if(fix(freq)~=freq), error('New frequency should be an integer'); end
  if(freq>p)
    error('The filter cutoff frequency should be smaller than current size');
  end
  if(nargin<3), useMat = true; end
  if(p>48), useMat = false; end

    %%-- Filtering
  FmatRet = [];
  if(useMat)
    if(size(matFreqs,2)>0)
      xIdx = find(matFreqs(1,:) == p);
      yIdx = find(matFreqs(2,xIdx) == freq);
    else
      xIdx = []; yIdx = [];
    end
    
    if(isempty(yIdx))
      Fmat{end+1} = getMat(p,freq);
      matFreqs(:,end+1) = [p;freq];
      xIdx = length(Fmat);
      yIdx = 1;
    end
    idx = xIdx(yIdx);
    ff = Fmat{idx}*f;
  else
    shc = shAna(f);
    shc((freq+1)^2+1:end,:) = 0;
    ff = shSyn(shc, false);
  end
  if( isreal(f)), ff = real(ff);end
end
  
function Fmat = getMat(p,freq)

  persistent queryUser

  np = 2*p*(p+1);
  printMsg('* Reading filtering matrix for p=%d to q=%d:',p,freq);
  fileName = ['Fmat' num2str(p) '_'];
  [Fmat, rFlag] = readData(fileName,freq,[np np]);
  if(~rFlag)
      printMsg(' Failed\n');
      printMsg('* Generating the filtering matrix for p=%d to q=%d.\n', ...
          p,freq);
      if(p>64)
          if(isempty(queryUser))
              queryUser = timeQuary([' Building the filter matrix for p=' num2str(p) ' may take ' ...
                  'a long time.\n Proceed with building? y/[n]?']);
              if(isempty(queryUser)), queryUser = 'n';end
          end
          if(strcmpi(queryUser,'n')), return; end
      end
    
      useMat = false; I = eye(np);
      Fmat  = filtersh(I, freq, useMat); writeData(fileName,freq, Fmat);
      printMsg('* Stored filter matrix for p=%d\n', p);
  else
    printMsg(' Successful\n');
  end
end