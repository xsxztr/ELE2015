function StokesMat = asmStokesP2P1(auxMesh,pde)
%% get mesh data
node=auxMesh.outnode;
elem=auxMesh.outelem;
elem2dof=auxMesh.outelem2dof;
edge=auxMesh.outedge;
% leftNode=auxMesh.leftNodeS;
% leftEdge=auxMesh.leftEdgeS;
% rightNode=auxMesh.rightNodeS;
% rightEdge=auxMesh.rightEdgeS;

bdFlagUz=auxMesh.bdFlagUz;
bdFlagUr=auxMesh.bdFlagUr;

% [node,elem,elem2dof,edge,elem2edge,leftNode,rightNode,leftEdge,rightEdge,bdFlagUz,bdFlagUr] = auxMesh{:};
N = size(node,1);  NT = size(elem,1);  Nu = N+size(edge,1);   Np = N;

%% Compute geometric quantities and gradient of local basis
% Dlambda=heateqMat.Dlambda;
% area = abs(simplexvolume(node,elem));
 [Dlambda,area] = gradbasis(node,elem);
 area = abs(area);
%% Assemble stiffness matrix for Laplace operator
% assemble the diffusion part \int r (\partial_ru_r\partial_rv_r+\partial_zu_r\partial_zv_r)
[lambda, w] = quadpts(3);
nQuad = size(lambda,1);
ii = zeros(21*NT,1); jj = zeros(21*NT,1); Aper = zeros(21*NT,1);
index = 0;
for i = 1:6
    for j = i:6
        Aij = 0;
        for p = 1:nQuad
            py = lambda(p,1)*node(elem(:,1),2) ...
                + lambda(p,2)*node(elem(:,2),2) ...
                + lambda(p,3)*node(elem(:,3),2);
            Aij = Aij + w(p)*dot(Dphi(p,i,lambda,Dlambda),Dphi(p,j,lambda,Dlambda),2).*py;
        end
        
        Aij = pde.eta*Aij;
        Aij = Aij.*area;
        ii(index+1:index+NT) = double(elem2dof(:,i));
        jj(index+1:index+NT) = double(elem2dof(:,j));
        Aper(index+1:index+NT) = Aij;
        index = index + NT;
    end
end

diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),Aper(diagIdx),Nu,Nu);
% A = spdiags(accumarray(ii(diagIdx),Aper(diagIdx),[Ndof 1]),0,Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),Aper(upperIdx),Nu,Nu);
A = A + AU + AU';
clear Aij ii jj Aper lambda w
%assimble \int u_rv_r/rdrdz
ii = zeros(21*NT,1); jj = zeros(21*NT,1); Aper = zeros(21*NT,1);
index = 0;
[lambda, w] = quadpts(4);
nQuad = size(lambda,1);
phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);
for i = 1:6
    for j = i:6
        Aij = 0;
        for p = 1:nQuad
            pr = lambda(p,1)*node(elem(:,1),2) ...
                + lambda(p,2)*node(elem(:,2),2) ...
                + lambda(p,3)*node(elem(:,3),2);
            Aij = Aij + w(p)*phi(p,i)*phi(p,j)./pr;
        end
        Aij = pde.eta*Aij;
        Aij = Aij.*area;
        ii(index+1:index+NT) = double(elem2dof(:,i));
        jj(index+1:index+NT) = double(elem2dof(:,j));
        Aper(index+1:index+NT) = Aij;
        index = index + NT;
    end
end
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
AR = sparse(ii(diagIdx),jj(diagIdx),Aper(diagIdx),Nu,Nu);
% A = spdiags(accumarray(ii(diagIdx),Aper(diagIdx),[Ndof 1]),0,Ndof,Ndof);
ARU = sparse(ii(upperIdx),jj(upperIdx),Aper(upperIdx),Nu,Nu);
AR=AR+ARU+ARU';
AR=AR+A;
A = blkdiag(A ,AR);
clear Aij ii jj Aper lambda w phi

%% Assemble the matrix for divergence operator
[lambda, w] = quadpts(3);
nQuad = size(lambda,1);
phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);
Dz = sparse(Np,Nu);
Dr = sparse(Np,Nu);
for i = 1:6
    for j = 1:3
        Dzij = 0;
        Drij = 0;
        
        for p = 1:nQuad
            pr = lambda(p,1)*node(elem(:,1),2) ...
                + lambda(p,2)*node(elem(:,2),2) ...
                + lambda(p,3)*node(elem(:,3),2);
            s = Dphi(p,i,lambda,Dlambda);
            Dzij = Dzij + w(p)*s(:,1).*lambda(p,j).*pr;
            Drij = Drij + w(p)*s(:,2).*lambda(p,j).*pr+w(p)* phi(p,i)*lambda(p,j);
            
        end
        Dz = Dz + sparse(elem(:,j),double(elem2dof(:,i)),Dzij.*area,Np,Nu);
        Dr = Dr + sparse(elem(:,j),double(elem2dof(:,i)),Drij.*area,Np,Nu);
    end
end
B = -[Dz Dr];
clear Dzij Drij Dz Dr

%% Fix peroidic BC
% pl=[leftNode(1:end-1);leftEdge; Nu+leftNode(2:end-1); Nu+leftEdge];
% pr=[rightNode(1:end-1);rightEdge; Nu+rightNode(2:end-1); Nu+rightEdge];
% A(:,pl)=A(:,pl)+A(:,pr);
% A(pl,:)=A(pl,:)+A(pr,:);
% B(:,pl)=B(:,pl)+B(:,pr);
[AD,BD,ufreeDofUr,ufreeDofUz,pDof] = getbdStokesP2P1;
StokesMat = struct('A',A,'B',B,'AD',AD,'BD',BD,'ufreeDofUr',ufreeDofUr,'ufreeDofUz',ufreeDofUz,'pDof',pDof);

%% ===============================================================
% subfunction getbdStokesP2P1
%% ===============================================================
    function [AD,BD,ufreeDofUr,ufreeDofUz,pDof] = getbdStokesP2P1
        
        %% Fix Dirichlet BC
        ufreeDofUr = (1:Nu)';
        ufreeDofUz = (1:Nu)';
        %% Part 1: Find Dirichlet dof and modify the matrix
        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedDofUr = false(Nu,1);
        isFixedDofUz = false(Nu,1);
        % case: bdFlag is not empty
        elem2edge = elem2dof(:,4:6)-N;
        isDirichletUr(elem2edge(bdFlagUr(:)==1)) = true;
        isFixedDofUr(edge(isDirichletUr,:)) = true;   % nodes of all D-edges
        isFixedDofUr(N + find(isDirichletUr')) = true;% dof on D-edges
        isDirichletUz(elem2edge(bdFlagUz(:)==1)) = true;
        isFixedDofUz(edge(isDirichletUz,:)) = true;   % nodes of all D-edges
        isFixedDofUz(N + find(isDirichletUz')) = true;% dof on D-edges
        %%
        fixedDofUr = find(isFixedDofUr);
        ufreeDofUr = find(~isFixedDofUr);
        fixedDofUz = find(isFixedDofUz);
        ufreeDofUz = find(~isFixedDofUz);
        % add peroidic nodes
        %fixedDofUr = [fixedDofUr; rightNode(2:end-1);rightEdge];
        %fixedDofUz = [fixedDofUz; rightNode(1:end-1);rightEdge];
        pDof = (1:Np-1)';
        
        %% Modify the matrix
        % Build Dirichlet boundary condition into the matrix AD by enforcing
        % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
        % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.
        bdidx = zeros(2*Nu,1);
        bdidx([fixedDofUz; Nu+fixedDofUr]) = 1;
        Tbd = spdiags(bdidx,0,2*Nu,2*Nu);
        T = spdiags(1-bdidx,0,2*Nu,2*Nu);
        AD = T*A*T + Tbd;
        BD = B*T;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction Dphi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function s = Dphi(p,i,lambda,Dlambda) % gradient of basis phi
        switch i
            case 1
                s = (4*lambda(p,1)-1).*Dlambda(:,:,1);
            case 2
                s = (4*lambda(p,2)-1).*Dlambda(:,:,2);
            case 3
                s = (4*lambda(p,3)-1).*Dlambda(:,:,3);
            case 4
                s = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
            case 5
                s = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
            case 6
                s = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end