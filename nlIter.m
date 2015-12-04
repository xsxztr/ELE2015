function [pos,neg,potent,u1,u2,press] = nlIter(potentOld,posOld,negOld,auxMesh,auxMesht,auxSphere,...
    pde,PoissMat,StokesMat,surfchar,svel,option)
%% Data Structure
% mesh
%node = auxMesh.node;
%elem = auxMesh.elem;
outnode = auxMesh.outnode;
outelem = auxMesh.outelem;
outnodeidx = auxMesh.outnodeidx;
%outelemidx = auxMesh.outelemidx;

outnodet = auxMesht.outnodet;
outelemt = auxMesht.outelemt;
%N = size(node,1);
%Nout = size(outnode,1);
Noutt = size(outnodet,1);
% pde

c = pde.c;
%% convex iteration%% remeber the last step iteration value for compare
 
err = 1;
 
%potentout=potentOld(outnodeidx);
while (err > option.nlItertol)
    %% solver poisson equation get potential
    
    f2 = zeros(Noutt,1); %% f2 is for the other ions in solvent 
    potentHalf = CylinderD(auxMesh,auxMesht,auxSphere,PoissMat,surfchar,posOld-negOld,f2,pde);
    potent = c*potentHalf+(1-c)*potentOld;
    %% compute the grad Psi
    
    potentout=potent(outnodeidx);
    V1 = gradu(outnode,outelem,potentout);
    
    net2=zeros(Noutt,1);
    V2=zeros(size(outelemt,1),2);
    [u1,u2,press] = StokesP2P1(auxMesh,auxMesht,StokesMat,pde,negOld-posOld,V1,net2,V2,svel);
    
    u=[u1,u2];
    poshalf =  getCharge(auxMesh,pde,potentout,u,1);
    neghalf = getCharge(auxMesh,pde,-potentout,u,2);
    
    pos=c*poshalf+(1-c)*posOld;
    neg=c*neghalf+(1-c)*negOld;
    
    
    
    %% compare the self-consistance
    erpot = max(abs(potent-potentOld));
    erpos = max(abs(pos-posOld));
    erneg = max(abs(neg-negOld));
    
    if (min(pos)<0||min(neg<0))
        
        fprintf('get negative distribution');
        break;
    end
    err = max([erpot,erpos,erneg ])/c;
    fprintf('nonlinear iteration err = %8.6g\n',err);
    posOld = pos;
    negOld = neg;
    potentOld = potent;
    
end
end