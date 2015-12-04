function PoissMat=asmPoiMat(auxMesh,pde)
%% mesh structure
node = auxMesh.node;
elem = auxMesh.elem;
innerelemidx= auxMesh.innerelemidx;
% leftNode = auxMesh.leftNode;
% rightNode = auxMesh.rightNode;
bdFlag = auxMesh.bdFlagPhi;

epss=pde.epss;
epsf=pde.epsf;
 N = size(node,1);  NT = size(elem,1);
 
 K=epsf*ones(NT,1);
 K(innerelemidx)=epss;
 
[Dlambda,area] = gradbasis(node,elem);
%% Assemble stiffness matrix
[lambda,w] = quadpts(2);
nQuad = size(lambda,1);

A = sparse(N,N);
for i = 1:3
    for j = i:3
        % $A_{ij}|_{\tau} = \int_{\tau}K\nabla \phi_i\cdot \nabla \phi_j dxdy$
        Aij = 0; 
        for p = 1:nQuad
            py = lambda(p,1)*node(elem(:,1),2) ...
                + lambda(p,2)*node(elem(:,2),2) ...
                + lambda(p,3)*node(elem(:,3),2);
        Aij = Aij+w(p)*(Dlambda(:,1,i).*Dlambda(:,1,j) + Dlambda(:,2,i).*Dlambda(:,2,j)).*py.*area;
        end
        Aij=K.*Aij;
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Aij; Aij],N,N);
        end
    end
end
clear K Aij lambda weight
%% Fix periodic boundary condition
% A(leftNode(1:end-1),:)=A(leftNode(1:end-1),:)+A(rightNode(1:end-1),:);
% A(:,leftNode(1:end-1))=A(:,leftNode(1:end-1))+A(:,rightNode(1:end-1));
[AD,freeNode] = getbd;
%% assimble 
PoissMat=struct('A',A,'AD',AD,'freeNode',freeNode,'area',area);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunction getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,freeNode] = getbd        
        %% Part 1: Find Dirichlet boundary nodes and modify the stiffness matrix
        % Find Dirichlet boundary nodes: fixedNode
        allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        Dirichlet = allEdge((bdFlag(:) == 1),:);
        isBdNode = false(N,1);
        isBdNode(Dirichlet(:)) = true;
        fixedNode = isBdNode;
        freeNode = find(~isBdNode);
         % add peroidic nodes
       % fixedNode = [fixedNode; rightNode(2:end-1)];
        bdidx = zeros(N,1);
        bdidx(fixedNode) = 1;
        Tbd = spdiags(bdidx,0,N,N);
        T = spdiags(1-bdidx,0,N,N);
        AD = T*A*T + Tbd;
end
end