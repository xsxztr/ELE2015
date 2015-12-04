function u=CylinderD(auxMesh,auxMesht,auxSphere,PoissMat,surfchar,f1,f2,pde)
%%% in the cylinder coordinate, here we take the length of cylinder as x,and
%%% the radius as y
%%%% p,n are in the \omega_f, but potential is in the whole domain 


%% mesh data
node = auxMesh.node;
%elem = auxMesh.elem;
outelem = auxMesh.outelem;
outnode = auxMesh.outnode;
outnodeidx = auxMesh.outnodeidx;

outelemt = auxMesht.outelemt;
outnodet = auxMesht.outnodet;
outnodeidxt = auxMesht.outnodeidxt;



% leftNode = auxMesh.leftNode;
% rightNode = auxMesh.rightNode;

C1 = 1/2;
A = PoissMat.A;
AD = PoissMat.AD;
ufreeNode = PoissMat.freeNode;
N = size(node,1);  
Nout=size(outnode,1); NTout = size(outelem,1);
Noutt=size(outnodet,1); NToutt = size(outelemt,1);
area =  abs(simplexvolume(outnode,outelem));
areat = abs(simplexvolume(outnodet,outelemt));
%% Assemble the right hand side
option.fquadorder = 3;   % default order
[lambda,weight] = quadpts(option.fquadorder);
phi = lambda;                 % linear bases
nQuad = size(lambda,1);
bt = zeros(NTout,3);
%%%% first part for H and OH%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    py = lambda(p,1)*outnode(outelem(:,1),2) ...
        + lambda(p,2)*outnode(outelem(:,2),2) ...
        + lambda(p,3)*outnode(outelem(:,3),2);
    fp = lambda(p,1)*f1(outelem(:,1)) ...
        + lambda(p,2)*f1(outelem(:,2)) ...
        + lambda(p,3)*f1(outelem(:,3));
    for i = 1:3
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i).*py.*fp*C1;
    end
end

clear py
bt = bt.*repmat(area,1,3);
bout = accumarray(outelem(:),bt(:),[Nout 1]);
b=zeros(N,1);
b(outnodeidx)=bout;


%%%%%% second part for salt %%%%%%%%%%%%%%%%%%%%%%%%%%
btt = zeros(NToutt,3);

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    py = lambda(p,1)*outnodet(outelemt(:,1),2) ...
        + lambda(p,2)*outnodet(outelemt(:,2),2) ...
        + lambda(p,3)*outnodet(outelemt(:,3),2);
    fpt = lambda(p,1)*f2(outelemt(:,1)) ...
        + lambda(p,2)*f2(outelemt(:,2)) ...
        + lambda(p,3)*f2(outelemt(:,3));
    for i = 1:3
        btt(:,i) = btt(:,i) + weight(p)*phi(p,i).*py.*fpt*C1;
    end
end

clear py
btt = btt.*repmat(areat,1,3);
boutt = accumarray(outelemt(:),btt(:),[Noutt 1]);
b(outnodeidxt)=b(outnodeidxt)+boutt;





%% Fix periodic boundary condition
%b(leftNode(1:end-1))=b(leftNode(1:end-1))+b(rightNode(1:end-1));
%% Set up boundary conditions
[b, u] = getbd(A,b);
%% Solve the system of linear equations
         tic; 
        u  =  AD\b;
%       err = norm(b -AD*u)/norm(b);
% 	time = toc;
% 	itStep=0;
% 	Ndof=size(AD,1);
% 	fprintf('Direct solver \n')
% 	fprintf('#dof: %8.0u,  #nnz: %8.0u,  iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
%     Ndof, nnz(A), itStep, err, time);

%% post-processes fix the periodic boundary condition
%u(rightNode(1:end-1))=u(leftNode(1:end-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [b,u] = getbd(A,b)
        u = zeros(N,1);
        
         Neumannout = auxSphere.outNeumann(2:end-1,:); 
         Neumann=outnodeidx(Neumannout);
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        option.gNquadorder = 2;
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),2);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),2) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),2);
            %gNp =surfchar*pde.phig_N(ppxy);
            for igN = 1:2
                %ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
                ge(:,igN) = ge(:,igN) + weightgN(pp)*surfchar*phigN(pp,igN)*ppxy;
            end
        end
        ge = ge.*repmat(el,1,2);
        b = b + accumarray(Neumann(:), ge(:),[N,1]); 
        
        
        
        
        
         isFixedNode = true(N,1);
        isFixedNode(ufreeNode)=false;
        fixedNode=find(isFixedNode);
        %% Part 2: Find boundary edges and modify the right hand side b
        u(fixedNode) = pde.phig_D(node(fixedNode,:));
        b = b - A*u;
        b(fixedNode) = u(fixedNode);
        %% ToDO: test Robin condition
    end % end of getbd
end