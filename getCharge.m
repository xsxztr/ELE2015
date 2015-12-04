function charge=getCharge(auxMesh,pde,phi,u,Iidx)
%% mesh stucture
%%%%%%%% new mesh structure %%%%%%%%%%5
node = auxMesh.outnode;
elem = auxMesh.outelem;

Elemcot = auxMesh.Elemcot;
%M = heateqMat.M;
N = size(node,1); NT = size(elem,1);

%%%%%%%% old mesh information %%%%%%%%%%%

%MOld = heateqMatOld.M;

alphae = 1;
bdFlag = auxMesh.bdFlagChar;
%  leftNode=auxMesh.leftNodeS;
% rightNode=auxMesh.rightNodeS;
% area = abs(simplexvolume(node,elem));
y =( node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2))/3;

if Iidx==1
    D=pde.Dp;
    char=pde.nH;
else
    D=pde.Dn;
    char=pde.nOH;
end

%% Assemble the EAFEM matrix
% [lambda,weight] = quadpts(2);
% nQuad = size(lambda,1);
A=sparse(N,N);
for i=1:3
    for j=1:3
        %         for p = 1:nQuad
        %             py= lambda(p,1)*node(elem(:,1),2) ...
        %                 + lambda(p,2)*node(elem(:,2),2) ...
        %                  + lambda(p,3)*node(elem(:,3),2);
        if i==j
            %use circle get the other two points of the element
            z=circle(i);
            %get vector \taue
            ep =  elem(:,z(1)) ;
            ip = elem(:,i);
            taue1=node(ip,:)-node(ep,:);
            %get the average of u on the edge E
            ue1=(u(ip(:),1).*taue1(:,1)+u(ip(:),2).*taue1(:,2)...
                +u(ep(:),1).*taue1(:,1)+u(ep(:),2).*taue1(:,2))./2;
            %get s for bernoulli function
            betatau1=(phi(ip)-phi(ep)-ue1)./alphae;
            
            ep = elem(:,z(2));
            ip =elem(:,i);
            taue2=node(ip,:)-node(ep,:);
            ue2=(u(ip(:),1).*taue2(:,1)+u(ip(:),2).*taue2(:,2)...
                +u(ep(:),1).*taue2(:,1)+u(ep(:),2).*taue2(:,2))./2;
            betatau2=(phi(ip)-phi(ep)-ue2)./alphae;
            
            Aii= y(:).*Elemcot(:,z(2))./2.*bornoli(-betatau1(:)).*alphae...
                + y(:).*Elemcot(:,z(1))./2.*bornoli(-betatau2(:)).*alphae;
            A=A+sparse(elem(:,i),elem(:,i),Aii,N,N);
        else
            ep =elem(:,i);
            ip =elem(:,j);
            taue=node(ip,:)-node(ep,:);
            % use findangle to get the theta_e
            k=findangle(i,j);
            ue=(u(ip(:),1).*taue(:,1)+u(ip(:),2).*taue(:,2)...
                +u(ep(:),1).*taue(:,1)+u(ep(:),2).*taue(:,2))./2;
            betatau=(phi(ip)-phi(ep)-ue)./alphae;
            
            Aji=- y(:).*Elemcot(:,k)./2.*bornoli(betatau(:)).*alphae;
            A=A+sparse(elem(:,j),elem(:,i), Aji,N,N);
            %         end
        end
    end
end
% clear lambda  weight py
%% Assemble right side term


%b= MOld*f;
b=zeros(N,1);


% option.fquadorder = 3;   % default order
% [lambda,weight] = quadpts(option.fquadorder);
% phi = lambda;                 % linear bases
% nQuad = size(lambda,1);
% bt = zeros(NT,3);
% for p = 1:nQuad
%     py= lambda(p,1)*node(elem(:,1),2) ...
%         + lambda(p,2)*node(elem(:,2),2) ...
%         + lambda(p,3)*node(elem(:,3),2);
%     fp = lambda(p,1)*f(elem(:,1),:) ...
%         + lambda(p,2)*f(elem(:,2),:) ...
%         + lambda(p,3)*f(elem(:,3),:);
%     for i = 1:3
%         bt(:,i) = bt(:,i) + weight(p).* py(:).*phi(p,i).*fp;
%     end
% end
% bt = bt.*repmat(area,1,3);
% b = accumarray(elem(:),bt(:),[N 1]);
A= A*D;
[charge,AD,b,freeNode]=getbd(A,b,bdFlag);
%% fix the periodic boundary condition
% A(:,leftNode)=A(:,leftNode)+A(:,rightNode);
% A(leftNode,:)=A(leftNode,:)+A(rightNode,:);
% b(leftNode)=b(leftNode)+b(rightNode);
% pflag=zeros(N,1);
% pflag(leftNode,1)=1;
% pflag(rightNode,1)=-1;
% npidx=find(pflag>-1);
% %% solve the equation directly
%  charge=zeros(N,1);
% tic;
% freeNode=npidx;

charge(freeNode) = AD(freeNode,freeNode)\b(freeNode);
% err = norm(b(freeNode)-A(freeNode,freeNode)*charge(freeNode))/norm(b(freeNode));
% time = toc;
% itStep=0;
% Ndof=size(freeNode,1);
% fprintf('Direct solver \n')
% fprintf('#dof: %8.0u,  #nnz: %8.0u,  iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
%     Ndof, nnz(A), itStep, max(err(end,:)), time);
%% post process
% charge(rightNode)=charge(leftNode);

%% deal with the dirichlet boundary condition
% %% boundary condition set
    function [charge,AD,b,freeNode]=getbd(A,b,bdFlag)
        charge = zeros(N,1);
        [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
        freeNode = find(~isBdNode);
        bdidx = zeros(N,1);
        bdidx(fixedNode) = 1;
        Tbd = spdiags(bdidx,0,N,N);
        T = spdiags(1-bdidx,0,N,N);
        AD = T*A*T + Tbd;
        
%         Neumann = bdEdge;
%         el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
%         
%         option.gNquadorder = 2;   % default order exact for linear gN
%         
%         [lambdagN,weightgN] = quadpts1(option.gNquadorder);
%         phigN = lambdagN;                 % linear bases
%         nQuadgN = size(lambdagN,1);
%         ge = zeros(size(Neumann,1),2);
%         for pp = 1:nQuadgN
%             % quadrature points in the x-y coordinate
%             ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
%                 + lambdagN(pp,2)*node(Neumann(:,2),:);
%             if Iidx==1
%                 gNp=pde.pg_N(ppxy);
%             else
%                 gNp=pde.ng_N(ppxy);
%             end
%             
%             for igN = 1:2
%                 ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
%             end
%         end
%         ge = ge.*repmat(el,1,2);
%         b = b + accumarray(Neumann(:), ge(:),[Ndof,1]);
        
         %% Part 2: Find boundary edges and modify the right hand side b
        charge(fixedNode) = char*pde.charg_D(node(fixedNode,:));
        b = b - A*charge;
        b(fixedNode) = charge(fixedNode);
    end
        
        
        
        
        %         fixedNode = [];
        %         allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        %         Dirichlet = allEdge((bdFlag(:) == 1),:);
        %         isBdNode = false(N,1);
        %         isBdNode(Dirichlet(:)) = true;
        %         fixedNode = find(isBdNode);
        %         freeNode = find(~isBdNode);
        %         bdidx = zeros(N,1);
        %         bdidx(fixedNode) = 1;
        %         Tbd = spdiags(bdidx,0,N,N);
        %         T = spdiags(1-bdidx,0,N,N);
        %         AD = T*A*T + Tbd;
        %         charge(fixedNode) = pde.Charg_D(node(fixedNode,:));
        %         b = b - A*charge;
        %         b(fixedNode) = charge(fixedNode);
        %     end
        function z=circle(point)
            switch point
                case 1
                    z=[2,3];
                case 2
                    z=[3,1];
                case 3
                    z=[1,2];
            end
        end
        function z=findangle(point1,point2)
            sumij=point1+point2;
            switch sumij
                case 3
                    z=3;
                case 4
                    z=2;
                case 5
                    z=1;
            end
        end
        function y=bornoli(s)
            
            y = s./(exp(s)-1);
            idx =  (s<10e-4);
            y(idx)=1-0.5*s(idx);
            idx =  (s<-50);
            y(idx)=-s(idx);
            idx =  (s>50);
            y(idx) = 0;
            idx =  (s==0);
            y(idx)=1;
        end
        
    end