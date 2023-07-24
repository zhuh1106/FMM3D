function [SMat, SpMat, DMat] = kernelDLap(S,type)
  persistent p Ywt A Wsph

  if ( nargin > 0 && ~isempty(S) )
    printMsg('  * Generating the direct double-layer matrix for p=%d.\n',S.p);
   
    np = length(S.cart.x);
    if (isempty(p) || S.p~=p)
      [Ywt A Wsph] = getPesistents(S.p);
      p = S.p;
    end
  
    W0 = S.geoProp.W./Wsph;
    DMat = zeros(np, np); SMat = DMat; SpMat = DMat;
    ind_ref = reshape((1:np),p+1,2*p);

    for k = 1:2*p
      ind = circshift(ind_ref,[0 k-1]);
      ind = ind(:);
      
      for j = 1:p+1
        %- Rotating the pole to (j,k)
        R = A{j}(ind,ind);
        xx = (R*S.cart.x)'; nx = (R*S.geoProp.nor.x)';
        yy = (R*S.cart.y)'; ny = (R*S.geoProp.nor.y)';
        zz = (R*S.cart.z)'; nz = (R*S.geoProp.nor.z)';
        
        W = ((R*W0).*Wsph)';
        ind_pole = (k-1)*(p+1) + j;

        %- The distance vector
        rx = S.cart.x(ind_pole) - xx;
        ry = S.cart.y(ind_pole) - yy;
        rz = S.cart.z(ind_pole) - zz;
        inv_rho = 1./sqrt(rx.^2 + ry.^2 + rz.^2);

        %- Normalizing
        rx = rx.*inv_rho;
        ry = ry.*inv_rho;
        rz = rz.*inv_rho;

        %- Laplace kernel and the scalar part
        g = (nx.*rx + ny.*ry + nz.*rz);
        g = g.*Ywt.*W.*inv_rho.*inv_rho;

        if nargin==1
        %- Single layer Matrix
        SMat(ind_pole,:) = (Ywt.*W.*inv_rho)*R;         
        %- Double layer Matrix
        DMat(ind_pole,:) = g*R;   
        %- Derivative of SLP
        SpMat(ind_pole,:) = -((S.geoProp.nor.x(ind_pole)*rx + S.geoProp.nor.y(ind_pole)*ry + S.geoProp.nor.z(ind_pole).*rz).*Ywt.*W.*inv_rho.*inv_rho)*R;
        else
            switch type
                case 'SMat'
                %- Single layer Matrix
                SMat(ind_pole,:) = (Ywt.*W.*inv_rho)*R;           
                case 'SpMat'
                %- Derivative of SLP
                SpMat(ind_pole,:) = -((S.geoProp.nor.x(ind_pole)*rx + S.geoProp.nor.y(ind_pole)*ry + S.geoProp.nor.z(ind_pole).*rz).*Ywt.*W.*inv_rho.*inv_rho)*R;    
                case 'DMat'
                %- Double layer Matrix
                DMat(ind_pole,:) = g*R;     
            end
        end
        
%         f_Temp = -((S.geoProp.nor.x(ind_pole)*rx + S.geoProp.nor.y(ind_pole)*ry + S.geoProp.nor.z(ind_pole).*rz).*W.*inv_rho.*inv_rho);
%         surf(reshape((g./Ywt)',p+1,2*p))
%         pause
      end
    end
    
  end

function [Ywt A Wsph] = getPesistents(p)

  printMsg('  * Generating the direct double-layer integration data for p=%d.\n',p);
  Ywt = (1/4/pi)*SingularWeights(p);
  u = gl_grid(p);
  Wsph = sin(u);
  np = 2*p*(p+1);

  A = cell(p+1,1);
  for idx = 1:p+1
    fname = ['RotMat' num2str(idx) '-'];
    A{idx} = readData(fname, p, [np np]);
    
    if(isempty(A{idx}))
      printMsg('  * Generating the %dth direct rotation matrices for p=%d.\n',idx,p-1);
      A{idx} = movePole(eye(np),u(idx),0);
      writeData(fname, p, A{idx});
    end
  end
  printMsg('  * Direct rotation matrices for p=%d were generated/read from file.\n',p);
end

end

function [fRot, sRot]= movePole(f,thetaTarg,lambdaTarg, isReal) 
% [cc shcNew]= movePole(f,thetaTarg,lambdaTarg) move the pole of the
%parameters to the target point with coordinate thetaTarg (elevation
%[-pi/2,pi/2]) and lambdaTarg (azimuth [0,2pi)).  MOVEPOLE returns the
%Cartesian coordinates of the points with now parametrization (fRot) and
%their corresponding spherical harmonic coefficients (shcRot).

  if(nargin==0), testMovePole(); return;end
  if(nargin<4), isReal = true; end
  [d1, d2] = size(f);
  N = (sqrt(2*d1+1)-1)/2;

  %% lambdaTarg rotation
  mMat = exp(1i*repmat(-N:N,N+1,1)*lambdaTarg);
  shc = repmat(shrinkShVec(mMat(:)),1,d2).*shAna(f);

  %% thetaTarg rotation
  %Parametrization of the equator in the rotated from. We need 2*N+2
  %equispaced points in [0,2pi).
  phiR = (0:2*N+1)*2*pi/(2*N+2) + pi/(2*N+2);
  %Finding the corresponding material point
  [phi, theta] = cart2sph(cos(thetaTarg)*cos(phiR),sin(phiR),-sin(thetaTarg)*cos(phiR));
  %Mapping phi to [0,2pi) 
  phi = mod(phi,2*pi)';
  theta = pi/2-theta'; 
  %Coefficients for G
  c1 = sin(thetaTarg)*cos(theta).*cos(phi) - cos(thetaTarg)*sin(theta);
  c2 = sin(thetaTarg)*sin(phi)./sin(theta);
  %Calculating f and g and finding their Fourier coefficients.
  for idx=1:d2
    for n=0:N    
      [Yn, Hn] = Ynm(n,[],theta(:),phi(:));
      coeff = repmat(shc(n^2+(1:2*n+1),idx).',2*N+2,1);
    
      Y2 = (c2*(1i*(-n:n))).*Yn;
      Hn = -repmat(c1,1,2*n+1).*Hn;
    
      f = sum(coeff.* Yn      ,2);
      g = sum(coeff.*(Hn + Y2),2);
    
      for m =-n:n
        ff(n+1,N+m+1) = 2*pi/(2*N+2)*sum(f.*exp(-1i*m*phiR'))/2/pi;
        gg(n+1,N+m+1) = 2*pi/(2*N+2)*sum(g.*exp(-1i*m*phiR'))/2/pi;
      end      
    end
  
    %Calculation new shc coeffiecents
    for n=0:N
      [P, Q] = Ynm(n,[],pi/2,0);
      for m=-n:n
        P0 = P(n+1+m);
        Q0 = Q(n+1+m);
        shcRot(n+1,N+m+1) = (ff(n+1,N+m+1)*P0 + gg(n+1,N+m+1)*Q0)/(P0^2+Q0^2);
      end
    end
    sRot(:,idx) = shcRot(:);
  end
  %% lambdaTarg rotation -- back 
  mMat = exp(-1i*repmat(-N:N,N+1,1)*lambdaTarg);
  sRot = repmat(mMat(:),1,d2).*sRot;
  sRot = shrinkShVec(sRot);
  
  %% Calculation point values 
  fRot = shSyn(sRot, isReal);
end

function val  = SingularWeights(n)

  LP_of_Y = zeros(1, n+1);
  for k = 0:n
    LP_of_Y(k+1) = 4*pi./(2*k+1)*Ynm(k,0,0,0);
  end

  [~, gwt]=g_grid(n+1); 
  wt = pi/n*repmat(gwt',1,2*n); 
  wt = wt(:)';
  [u, v] = gl_grid(n);
  Yf = 0; 
  for j = 0:n
    Yf = Yf + LP_of_Y(j+1)*(reshape(Ynm(j,0,u,v),2*n*(n+1),1)');
  end
  val = wt.*Yf./cos(u(:)'/2);
end

