function MMC3D_2(DL,DW,DH,nelx,nely,nelz,xInt,yInt,zInt,vInt,volfrac)
% DL=64 ;
% DW=8 ;
% DH=32 ;
% 
% nelx=128;
% nely=16;
% nelz=64 ;
% 
% xInt=8 ;
% yInt=4 ;
% zInt=8 ;
% 
% vInt=[12,2.5,2.0,0,atan(1),0] ;
% volfrac=0.3;
E = 1.0;	
nu = 0.3;	 % Young's modulus, Poisson ratio
dgt0 = 5;            
scl = 1.0;	  % significant digit of sens., scale factor for obj.
p = 6;	
lmd = 100;	 % power of super ellipsoid , power of KS aggregation
iter = 1;            
maxiter = 500;                                                                                          % initial and maximum number of iterations 
objVr5 = 1.0;

   %-------------------------SEC 2):   SETTING OF FE DISTRIBUTION
nEle = nelx*nely*nelz ;	%  number of finite elements
nNod = (nelx+1)*(nely+1)*(nelz+1) ;	
nNodfc = (nelx+1)* (nely+1) ;     %  numbers of all nodes and nodes in each layer
nDof = 3*nNod ;	 % number of degree of freedoms
EL = DL/nelx ;            
EW = DW/nely ;         
EH = DH/nelz ;	% length. width. height of each finite element
minSz = min([EL, EW, EH]) ;	% minimum size of finite elements
alpha=1e-3 ;     
epsilon = 0.2;	% void density

[Ke] = Ke_tril (E, nu, EL, EW, EH); 		 % non-zero upper triangular of ele. stiffness
KE (tril (ones(24))==1 ) = Ke' ;

KE = reshape (KE, 24, 24) ;	
KE = KE + KE' - diag(diag(KE)) ;

nodMat = int32(reshape(1:nNod, 1+nelx, 1+nely, 1+nelz));	% full elemental stiffness matrix
edofVec =reshape(3*nodMat(1:nelx, 1:nely, 1:nelz) , nEle, 1);	% matrix of nodes number  (int32)
edofMat = edofVec + int32([3*nNodfc+[[3*nelx+[4 5 6 1 2 3]] -2 -1 0 1 2 3] 3*nelx+[4 5 6 1 2 3] -2 -1 0 1 2 3]);	% connectivity matrix
eleNodesID = edofMat(: ,[3:3:24])./3;
[sI, sII] = deal ([]);
for j = 1:24
    sI =cat(2, sI, j : 24);
    sII = cat(2, sII, repmat(j, 1, 24-j+1));
end

[iK, jK] = deal(edofMat(:, sI)', edofMat (:, sII)');

Iar0 =sort([iK(:),jK(:)], 2, 'descend') ; 
clear iK jK	          % reduced assembly indexing

[x, y, z] = meshgrid(EL*[0:nelx], EW*[0:nely], EH*[0:nelz]);	         % coordinates of nodal points
LSgrid.x = permute(x, [2,1,3]) ; LSgrid.y = permute(y, [2,1,3]) ; LSgrid.z = permute(z, [2,1,3]) ;
volNod = sparse(double(eleNodesID(:)), 1, 1/8);

   %---------------------------SEC 3):   LOADS. DISPLACEMENT BOUNDARY CONDITIONS (3D cantilever beam example)
[jN,kN] = meshgrid(1:nely+1,1:nelz+1); 
fixNd1 = [(kN-1)*nNodfc+(jN-1)*(nelx+1)+1];
[iN,kN] = meshgrid(1:nelx+1,1:nelz+1); 
fixNd2 = [kN*nNodfc-iN+1];
fixNd3 = nelx + 1;
fixDof = [3*fixNd1(:)-2; 3*fixNd2(:)-1; 3*fixNd3(:);
3*fixNd3(:)-1; 3*fixNd3(:)-2];
fixEle = [nelx];
freeDof = setdiff([1:nDof],fixDof);
loadNd = nNod - nelx; loadDof = 3*loadNd;
loadEle = [nEle-nelx+1];
F = fsparse(loadDof',1,-1/4,[nDof,1]);

   %---------------------------------------------- SEC 4): INITIAL SETTING OF COMPONENTS

x0 = xInt : 2*xInt:DL; 
y0 = yInt: 2*yInt : DW; 
z0= zInt : 2*zInt:DH; % coordinates of inital components centre

xn = length(x0) ;  
yn = length(y0);  
zn=length(z0) ;% num. pf components along x,y,z

% x0 = kron(x0,ones(1,2*yn*zn));                                         % full xo vector
% y0 = repmat(repmat(y0,1,zn), 1, 2*xn ) ;                               % full y0 vector
% z0 = repmat(reshape(repmat(z0,2,1),1,4,[]),1,xn);                        % full zo vector

x0 = kron(x0,ones(1,4*yn*zn));
y0 = repmat(repmat(y0,1,2*zn),1,2*xn);
z0 = repmat(reshape(repmat(z0,4,1),1,8),1,xn);

N= length(x0); 
l1=repmat(vInt(1),1,N) ;
l2=repmat(vInt(2),1,N) ;
l3=repmat(vInt(3),1,N) ;
                                          % vector of half height
alp = repmat(vInt (4), 1, N);                                             % vector of alpha
bet = repmat([1 -1 1 -1 1 -1 1 -1]*vInt(5),1, N/8);                    % vector of beta
gam = repmat([1 1 -1 -1 1 1 -1 -1]*vInt(6),1, N/8);                     % vector of gamma
dd = [x0 ; y0 ; z0 ; l1 ; l2 ; l3 ; alp ; bet ; gam] ;                     % design of variable vector
nDsvb = length(dd(:));                                                         % number of all design variables
nEhcp = nDsvb/N;                                                                 %number of design variables each components
actComp= [1 : N];                                                              % initial set of active components
actDsvb= [ 1 : nDsvb];                                                         % initial set of active design variables
nNd = 0;    PhiNd = [];                                                       % number of non design patches and its TDF  matrix
allPhi = [zeros(nNod,N) PhiNd];                                                  % initialized TDF matrix




% ---------------------------------------------SEC 5): SETTING OF MMA
m = 1; c = 1000*ones( m,1);  d = zeros(m ,1 ) ;                              % MMA parameters 
a0=1; a= zeros(m,1);
xval = dd(:);     xold1 = xval;    xold2 = xval;
xmin = [ 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; -pi ; -pi ; -pi ] ;             %lower bounds
xmax = [ DL; DW ; DH ; sqrt( DL^2+ DW^2+ DH^2)/ 2*[1; 1; 1]; pi; pi; pi];      % upper bounds
xmin = repmat(xmin,N,1);     xmax = repmat (xmax,N,1);
low = xmin;   upp = xmax;

% ---------------------- SEC 6) : OPTIMIZATION LOOP
while objVr5>1e-4 && iter<=maxiter
    %-----------------------------LP 1): Generating TDFs and their
    %derivatives
    allPhiDrv = sparse(nNod, nDsvb);
    for i = actComp
        [allPhi, allPhiDrv, xval, actComp, actDsvb] = ...
            calc_Phi(allPhi, allPhiDrv, xval, i,LSgrid,p,nEhcp,epsilon,actComp,actDsvb,minSz);
    end
    allPhiAct = [allPhi(:,actComp) PhiNd] ;
    temp = exp(lmd*allPhiAct);
    Phimax = max(-1e3, log(sum(temp,2))/lmd);
    allPhiDrvAct = allPhiDrv(:,actDsvb);
    Phimaxdphi = kron(temp(:,1:length(actComp))./(sum(temp,2)+eps),ones(1,nEhcp));
    PhimaxDrvAct = Phimaxdphi.*allPhiDrvAct;
    %--------------------------------------LP 2): Plotting current design
     clf; h = patch(isosurface(x,y,z,permute(reshape(Phimax,nelx+1,nely+1,nelz+1),[2,1,3]),0));
     h1 = patch(isocaps(x,y,z,permute(reshape(Phimax,nelx+1,nely+1,nelz+1),[2,1,3]),0));
     set(h,'FaceColor','red','EdgeColor','none','facealpha',1);set(h1,'FaceColor','interp','EdgeColor','none');
     colormap([1 0 0]); isonormals(x,y,z,permute(reshape(Phimax,nelx+1,nely+1,nelz+1),[2,1,3]),h); lighting phong;
     view(3); axis image; axis([0,DL,0,DW,0,DH]); camlight right; pause(1e-1);
    %--------------------------------------LP 3) : Finite Element Analysis
    H = Heaviside(Phimax,alpha,epsilon);
    den = sum(H(eleNodesID),2)/8;
    U = zeros(nDof,1);
    nAct = length(actComp) + nNd;
    [strct, loadPth] = srch_ldpth(nAct, allPhiAct,Phimax,epsilon,eleNodesID, loadEle, fixEle);
    if strct == 1
        if length(loadPth) == nAct 
            denSld = den;
        else
            PhimaxSld = max(-1e3, log(sum(exp(lmd*allPhiAct(:,loadPth)),2))/lmd);
            HSld = Heaviside(PhimaxSld, alpha, epsilon);
            denSld = sum(HSls(eleNodesID),2)/8;
        end
        eleLft = setdiff([1:length(denSld)],find(denSld<alpha+eps));
        edofMatLft = edofMat(eleLft, : );
        freedofLft = setdiff(edofMatLft, fixDof);
        [iK1, jK1] = deal(edofMatLft(:,sI)',edofMatLft(:,sII)');
        Iar = sort([iK1(:),jK1(:)],2,'descend');clear iK1 jK1;
        sK = reshape(Ke(:)*denSld(eleLft)', length(Ke)*length(eleLft),1);
        K = fsparse(Iar(:,1), Iar(:,2),sK,[nDof, nDof]); K=K+K'-diag(diag(K));
        K= K+ fsparse([1:nDof],[1:nDof],eps*ones(1,nDof),[nDof,nDof]);
        U(freedofLft) = K(freedofLft, freedofLft)\F(freedofLft);
    else
        disp('WARNING!!! NO load transmission path is founded!!!');
        sK = reshape(Ke(:)*den(:)',length(Ke)*nEle,1);
        K = fsparse(Iar0(:,1),Iar0(:,2),sK,[nDof,nDof]); K=K+K'-diag(diag(K));
        U(freeDof) = K(freeDof,freeDof)\F(freeDof);
    end
    f0val = F'*U/scl; OBJ(iter)=f0val*scl;
    fval = sum(den)*EL*EW*EH/(DL*DW*DH) - volfrac; CONS(iter) = fval + volfrac;
    
    %--------------------------------- LP 4): Sensitivity Analysis
    df0dx = zeros(1,nDsvb); dfdx = zeros(1,nDsvb);
    delta_H = 3*(1-alpha)/(4*epsilon)*(1-Phimax.^2/(epsilon^2));
    delta_H(abs(Phimax)>epsilon) = 0;
    energy = sum((U(edofMat)*KE).*U(edofMat),2);
    sEner = energy*ones(1,8)/8;
    engyNod = sparse(double(eleNodesID(:)),1,sEner(:));
    %df0dx(actDsvb) = -(engyNod.*delta_H)'*PhimaxDrvAct;
    dfdx(actDsvb) = (volNod.*delta_H)'*PhimaxDrvAct*EL*EW*EH/(DL*DW*DH);
    dgt = dgt0 - floor(log10([max(abs(df0dx(:))) max(abs(dfdx(:)))]));
    df0dx = round(df0dx*10^dgt(1))/10^dgt(1)/scl;
    dfdx = round(dfdx*10^dgt(2))/10^dgt(2);
    %--------------------LP 5): Updating design variables
     [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m, nDsvb, iter, xval(:),xmin,xmax,xold1,xold2, f0val, df0dx', fval, dfdx, low,upp,a0,a,c,d);
     xold2 = xold1;   xold1 = xval(:);   xval = xmma;                % design variable's update
     if iter>=5 && fval/volfrac<11-4
         objVr5 = abs(max(abs(OBJ(iter-4:iter)-mean(OBJ(iter-4:iter))))/mean(OBJ(iter-4:iter)));
     end
     disp([' It.: ' sprintf('%4i\t',iter) ' Obj.: ' sprintf('%6.3f\t',f0val*scl) ' Vol.: ' sprintf('%6.4f\t',fval) 'ch.:' sprintf('%6.4f\t',objVr5)]);
     iter = iter + 1;
end



