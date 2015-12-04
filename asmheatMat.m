function heateqMat=asmheatMat(auxMesh,auxMesht)
%% mesh data structure
node=auxMesh.outnode;
elem=auxMesh.outelem;
%N=size(node,1);  
 
 
%% Compute geometric quantities and gradient of local basis
 [Dlambda,area] = gradbasis(node,elem);
%  
% %% Assemble the mass matrix
% M=sparse(N,N);
% [lambda,weight] = quadpts(3);
% nQuad = size(lambda,1);
% for i = 1:3
%     for j = 1:3
%         Sij=0;
%         for p = 1:nQuad
%             py  = lambda(p,1)*node(elem(:,1),2) ...
%                 + lambda(p,2)*node(elem(:,2),2) ...
%                 + lambda(p,3)*node(elem(:,3),2);
%             Sij=Sij+weight(p)*lambda(p,i)*lambda(p,j).*py.*area;
%         end
%         M=M+sparse(elem(:,i),elem(:,j),Sij,N,N);
%     end
% end




%% mesh out of the trap %%%%%%%%%%%%%%%

 nodet=auxMesht.outnodet;
 elemt=auxMesht.outelemt;
% Nt=size(nodet,1);  
%  
%  
% %% Compute geometric quantities and gradient of local basis
[Dlambdat,areat] = gradbasis(nodet,elemt);
%  
% %% Assemble the mass matrix
% Mt=sparse(Nt,Nt);
% [lambda,weight] = quadpts(3);
% nQuad = size(lambda,1);
% for i = 1:3
%     for j = 1:3
%         Sij=0;
%         for p = 1:nQuad
%             py  = lambda(p,1)*nodet(elemt(:,1),2) ...
%                 + lambda(p,2)*nodet(elemt(:,2),2) ...
%                 + lambda(p,3)*nodet(elemt(:,3),2);
%             Sij=Sij+weight(p)*lambda(p,i)*lambda(p,j).*py.*areat;
%         end
%         Mt=Mt+sparse(elemt(:,i),elemt(:,j),Sij,Nt,Nt);
%     end
% end
% 





 
clear Sij lambda weight
heateqMat=struct('Dlambda',Dlambda,'area',area,'Dlambdat',Dlambdat,'areat',areat);
end