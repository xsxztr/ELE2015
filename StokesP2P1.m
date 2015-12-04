function [u1,u2,p] = StokesP2P1(auxMesh,auxMesht,StokesMat,pde,net1,V1,net2,V2,svel)
%
%% Construct Data Structure
nodeall=auxMesh.node;
elemall=auxMesh.elem;
node=auxMesh.outnode;
elem=auxMesh.outelem;
nodet=auxMesht.outnodet;
elemt=auxMesht.outelemt;
outelemidxt=auxMesht.outelemidxt;
outelemidx=auxMesh.outelemidx;
elem2dof=auxMesh.outelem2dof;
edge=auxMesh.outedge;
 
 



 
A = StokesMat.A;
B = StokesMat.B;
AD = StokesMat.AD;
BD = StokesMat.BD;
ufreeDofUr = StokesMat.ufreeDofUr;
ufreeDofUz = StokesMat.ufreeDofUz;
pDof = StokesMat.pDof;

N = size(node,1);  NT = size(elem,1);  Nu = N+size(edge,1);   Np = N;
Nt = size(nodet,1);  NTt = size(elemt,1);  
Nall = size(nodeall,1);  NTall = size(elemall,1); 


area = abs(simplexvolume(node,elem));
areat =abs(simplexvolume(nodet,elemt));
%% Compute geometric quantities and gradient of local basis
% quadrature points in the barycentric coordinate
[lambda,weight] = quadpts(4);
nQuad = size(lambda,1);

%%   % basis values at quadrature points
phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);
%% Assemble right hand side
%%%%% fist part for H and OH %%%%%%%%%%%%%%5
ft1 = zeros(NT,6);
ft2 = zeros(NT,6);
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    netxy = lambda(p,1)*net1(elem(:,1),:) ...
        + lambda(p,2)*net1(elem(:,2),:) ...
        + lambda(p,3)*net1(elem(:,3),:);
    py= lambda(p,1)*node(elem(:,1),2) ...
        + lambda(p,2)*node(elem(:,2),2) ...
        + lambda(p,3)*node(elem(:,3),2);
    % function values at quadrature points
    fp = repmat(netxy,1,2).*V1;
    % evaluate fp outside.
    for j = 1:6
        ft1(:,j) = ft1(:,j) + fp(:,1).*phi(p,j)*weight(p).*py(:);
        ft2(:,j) = ft2(:,j) + fp(:,2).*phi(p,j)*weight(p).*py(:);
    end
end
ft1 = ft1.*repmat(area,1,6);
ft2 = ft2.*repmat(area,1,6);
f1 = accumarray(elem2dof(:),ft1(:),[Nu 1]);
f2 = accumarray(elem2dof(:),ft2(:),[Nu 1]);

clear py ft1 ft2 fp

%%%%%%second part for salt %%%%%%%%%%%%%%%%%
ft1 = zeros(NTt,6);
ft2 = zeros(NTt,6);
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    netxyt = lambda(p,1)*net2(elemt(:,1),:) ...
        + lambda(p,2)*net2(elemt(:,2),:) ...
        + lambda(p,3)*net2(elemt(:,3),:);
    py= lambda(p,1)*nodet(elemt(:,1),2) ...
        + lambda(p,2)*nodet(elemt(:,2),2) ...
        + lambda(p,3)*nodet(elemt(:,3),2);
    % function values at quadrature points
    fp = repmat(netxyt,1,2).*V2;
    % evaluate fp outside.
    for j = 1:6
        ft1(:,j) = ft1(:,j) + fp(:,1).*phi(p,j)*weight(p).*py(:);
        ft2(:,j) = ft2(:,j) + fp(:,2).*phi(p,j)*weight(p).*py(:);
    end
end
ft1 = ft1.*repmat(areat,1,6);
ft2 = ft2.*repmat(areat,1,6);

ft1all=zeros(NTall,6);
ft2all=zeros(NTall,6);
ft1all(outelemidxt,:)=ft1;
ft2all(outelemidxt,:)=ft2;

 
ft1t=ft1all(outelemidx,:);
ft2t=ft2all(outelemidx,:);


f1 = f1+accumarray(elem2dof(:),ft1t(:),[Nu 1]);
f2 = f2+accumarray(elem2dof(:),ft2t(:),[Nu 1]);



clear py ft1 ft2 fp ft1all ft2all ft1t ft2t






%% Fix Periodic BC
% perLUz = [leftNode(1:end-1);leftEdge];
% perRUz = [rightNode(1:end-1);rightEdge];
% perLUr = [leftNode(2:end-1);leftEdge];
% perRUr = [rightNode(2:end-1);rightEdge];
% f1(perLUz)=f1(perLUz)+f1(perRUz);
% f2(perLUr)=f2(perLUr)+f2(perRUr);
% f1(perRUz)=0;
% f2(perRUr)=0;

%%
[f,g,u,p] = getbdStokesP2P1;
% pflag=zeros(2*Nu,1);
% pflag(perL)=1;
% pflag(Nu+perL)=1;
% pflag(perR)=-1;
% pflag(perR+Nu)=-1;
% npidx=find(pflag>-1);

%% Solver
%[u(npidx),p] = adgs(AD(npidx,npidx),BD(:,npidx),f(npidx),g,1e3*N);
%[u,p] = adgs(AD,BD,f,g,1e3*N);
bigA = [AD BD'; BD sparse(Np,Np)];
bigF = [f; g];
bigu = [u;p];
bigF = bigF - bigA*bigu;
bigFreeDof = [ufreeDofUz; Nu+ufreeDofUr; 2*Nu+pDof];
bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
u = bigu(1:2*Nu);
p = bigu(2*Nu+1:end);
u1 = u(1:N);
u2 = u(Nu+1:Nu+N);

%% Post-process
if length(pDof)~=Np % p is unique up to a constant
    % impose the condition int(r*p)=0
    r=node(:,2);
    c = sum(mean(r(elem).*p(elem),2).*area)/sum(area);
    p = p - c;
end

% u(perRUz)=u(perLUz);
% u(perRUr+Nu)=u(perLUr+Nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdStokesP2P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [f,g,u,p] = getbdStokesP2P1
        %% Boundary condition of Stokes equation: P2-P0 elements
        
        %% Initial set up
        g = zeros(Np,1);
        u = zeros(2*Nu,1);
        p = zeros(Np,1);
        f = [f1; f2];
        isFixedDofUr = true(Nu,1);
        isFixedDofUr(ufreeDofUr)=false;
        fixedDofUr=find(isFixedDofUr);
         isFixedDofUz = true(Nu,1);
        isFixedDofUz(ufreeDofUz)=false;
        fixedDofUz=find(isFixedDofUz);
        % Dirichlet boundary conditions
      
            ur = zeros(Nu,1);
            uz = zeros(Nu,1);
            idxUr = (fixedDofUr > N);              % index of edge dof
            idxUz = (fixedDofUz > N);    
            uDUr = pde.Urg_D(node(fixedDofUr(~idxUr),:));  % bd value at vertex dofs
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uDUz =svel* pde.Uzg_D(node(fixedDofUz(~idxUz),:));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ur(fixedDofUr(~idxUr)) = uDUr;
            uz(fixedDofUz(~idxUz)) = uDUz;
            bdEdgeIdxUr = fixedDofUr(idxUr)-N;
            bdEdgeIdxUz = fixedDofUz(idxUz)-N;
            bdEdgeMidUr = (node(edge(bdEdgeIdxUr,1),:)+node(edge(bdEdgeIdxUr,2),:))/2;
            bdEdgeMidUz = (node(edge(bdEdgeIdxUz,1),:)+node(edge(bdEdgeIdxUz,2),:))/2;
            uDUr = pde.Urg_D(bdEdgeMidUr);         % bd values at middle points of edges
            %%% svel is the velocity of sphere, the velocity of fluid
            %%% should be  equal to svel.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
            uDUz = svel*pde.Uzg_D(bdEdgeMidUz); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ur(fixedDofUr(idxUr)) = uDUr;
            uz(fixedDofUz(idxUz)) = uDUz;
            u = [uz; ur]; % Dirichlet bd condition is built into u
            f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
            g = g - B*u;  % the right hand side
            f(fixedDofUz) = uz(fixedDofUz);
            f(fixedDofUr+Nu) = ur(fixedDofUr);
       
    end % end of function getbdStokesP2P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
