function u=Cylinder(auxMesh,auxMesht,auxSphere,heateqMat,PoissMat,surfchar,f1,f2,pde)
%%% in the cylinder coordinate, here we take the length of cylinder as x,and
%%% the radius as y
%%%% p,n are in the \omega_f, but potential is in the whole domain


%% mesh data
node = auxMesh.node;
elem = auxMesh.elem;
N = size(node,1);
outnodeidx = auxMesh.outnodeidx;
outelem = auxMesh.outelem;
outnode = auxMesh.outnode;


outelemt = auxMesht.outelemt;
outnodet = auxMesht.outnodet;
outnodeidxt = auxMesht.outnodeidxt;
C1 = 1/2;
A = PoissMat.A;
bdFlag = auxMesh.bdFlagPhi;
% AD = PoissMat.AD;
% freeNode = PoissMat.freeNode;

Nout=size(outnode,1); NTout = size(outelem,1);
Noutt=size(outnodet,1); NToutt = size(outelemt,1);
area = heateqMat.area;
areat = heateqMat.areat;
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
% b=zeros(Nout,1);
% b=bb(outnodeidx); 
 
%% Set up boundary conditions
[AD,b,u,freeNode] = getbd(A,b);
%% Solve the system of linear equations 
u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 function [AD,b,u,freeNode] = getbd(A,b)
        Ndof=N;
        u=zeros(Ndof,1);
        fixedNode = []; freeNode = []; 
        [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
        freeNode = find(~isBdNode);
        bdidx = zeros(Ndof,1); 
        bdidx(fixedNode) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
        
        Neumann = auxSphere.outNeumann; 
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        option.gNquadorder = 2;
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),2);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            %gNp =surfchar*pde.phig_N(ppxy);
            for igN = 1:2
                %ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
                ge(:,igN) = ge(:,igN) + weightgN(pp)*surfchar*phigN(pp,igN);
            end
        end
        ge = ge.*repmat(el,1,2);
        b = b + accumarray(Neumann(:), ge(:),[Ndof,1]); 
        
        u(fixedNode) = pde.phig_D(node(fixedNode,:));
        b = b - A*u;
 end
% 
%     function [b,u] = getbd(A,b)
%         u = zeros(Nout,1);
%         
%         fixedNode=PoissMat.fixedNode;
%         Neumann=PoissMat.bdEdge;
%         %% Neumann boundary condition
%         
%         el = sqrt(sum((outnode(Neumann(:,1),:) - outnode(Neumann(:,2),:)).^2,2));
%         
%         option.gNquadorder = 2;   % default order exact for linear gN
%         
%         [lambdagN,weightgN] = quadpts1(option.gNquadorder);
%         phigN = lambdagN;                 % linear bases
%         nQuadgN = size(lambdagN,1);
%         ge = zeros(size(Neumann,1),2);
%         for pp = 1:nQuadgN
%             % quadrature points in the x-y coordinate
%             ppxy = lambdagN(pp,1)*outnode(Neumann(:,1),:) ...
%                 + lambdagN(pp,2)*outnode(Neumann(:,2),:);
%             gNp =surfchar*pde.phig_N(ppxy);
%             for igN = 1:2
%                 ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
%             end
%         end
%         ge = ge.*repmat(el,1,2);
%         
%         b = b + accumarray(Neumann(:), ge(:),[Nout,1]);
%          %showsolution(outnode,outelem,b);
%         
%         
%         
%         
%         %% Part 2: Find boundary edges and modify the right hand side b
%         u(fixedNode) = pde.phig_D(outnode(fixedNode,:));
%        % showsolution(outnode,outelem,u);
%         b = b - A*u;
%         b(fixedNode) = u(fixedNode);
%         %% ToDO: test Robin condition
%     end % end of getbd
end