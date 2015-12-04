clear all
clc
%% parameters
%%%%%%%%% these parameters can not be changed for the water system%%%%%%%%%%
par.dy = 0.003;%% the width of trap
par.eta = 24.7787; %% viscosity of fluid

par.L = 30;%% length of domain[-L, L]
par.epss=2.4;%% dielectric constant of sphere
par.epsf=80;%% dielectric constant of fluid
par.R = 2; %% radius of sphere
par.R0 = 30;%% radius of channel
par.phi0=0.0;%% hight of trap
par.theta=pi/10000; %% step size of the sphere surface  mesh


par.C1 = 0.5668; %% dimesion number before \frac{\partial ^2A}{\partal t^2}
par.C2 = 0;%% dimesion number before k
SpotOld=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%these parameters can be changed for different ions %%%%%%%%%%%%%%

par.Vmax=0.3;%% potential drop added on the both side  of channel
%%%%diffusion coefficient is dimensionless by the diffution coefficient of
%%%%H;
par.Dn=1;%% diffusion coefficient of OH
par.Dp=1;%% diffusion coefficient of H
par.Dk=0.2103;%%diffusion coefficient of  salt K
par.Dcl=0.2178;%% diffusion coefficient of cl

%%%% charge density is dimensionless by the 950.  If adding salt, nk and ncl should not be zeros.
%%%% if adding acid, the nH is equal to the density of H in the water and acid, and nk is zero, ncl
%%% % is equal to the density of negative ion of acid; If adding Alkali, the nOH  is equal to the density of  oH in the water and Alkali and ncl is zero, nk
%%% % is equal to the density of positive ion of aciddetails can be found in the note.

par.nH=1;%% density of H;
par.nOH=1;%% density of OH
par.nK=1; %% density of positive ions
par.ncl=1; %% density of negative ions
%par.surchar=1;
%% mesh data
load('auxMesh.mat');
load('auxMesht.mat');
node = auxMesh.node;
elem = auxMesh.elem;
outnode=auxMesh.outnode;
outelem=auxMesh.outelem;
innernode=auxMesh.innernode;
innerelem=auxMesh.innerelem;
outnodeidx=auxMesh.outnodeidx;
outelemidx=auxMesh.outelemidx;
innernodeidx=auxMesh.innernodeidx;
innerelemidx=auxMesh.innerelemidx;


outnodet=auxMesht.outnodet;
outelemt=auxMesht.outelemt;
outnodeidxt=auxMesht.outnodeidxt;
outelemidxt=auxMesht.outelemidxt;


% %% mesh  revise
% idx=(abs(node(:,2))<0.0001);
% node(idx,2)=0;
%
% idx=(abs(outnode(:,2))<0.0001);
% outnode(idx,2)=0;
%
% idx=(abs(outnodet(:,2))<0.0001);
% outnodet(idx,2)=0;
%
% idx=(abs(innernode(:,2))<0.0001);
% innernode(idx,2)=0;


%% information for EAFEM in the $\Omega_s$ area
area = simplexvolume(outnode,outelem);
NTout=size(outelem,1);
%edge vector of each element 23, 31,12,
vector(:,1:2)=outnode(outelem(:,3),1:2)-outnode(outelem(:,2),1:2);
vector(:,3:4)=outnode(outelem(:,1),1:2)-outnode(outelem(:,3),1:2);
vector(:,5:6)=outnode(outelem(:,2),1:2)-outnode(outelem(:,1),1:2);
%elem cot
Elemcot =zeros(NTout,3);
Elemcot(:,1)=-(vector(:,3).*vector(:,5)+vector(:,4).*vector(:,6))./(2*area(:));
Elemcot(:,2)=-(vector(:,1).*vector(:,5)+vector(:,2).*vector(:,6))./(2*area(:));
Elemcot(:,3)=-(vector(:,3).*vector(:,1)+vector(:,4).*vector(:,2))./(2*area(:));



%% information for EAFEM in the $\Omega_st$ area
areat = abs(simplexvolume(outnodet,outelemt));
NToutt=size(outelemt,1);
%edge vector of each element 23, 31,12,
vectort(:,1:2)=outnodet(outelemt(:,3),1:2)-outnodet(outelemt(:,2),1:2);
vectort(:,3:4)=outnodet(outelemt(:,1),1:2)-outnodet(outelemt(:,3),1:2);
vectort(:,5:6)=outnodet(outelemt(:,2),1:2)-outnodet(outelemt(:,1),1:2);
%elem cot
Elemcott =zeros(NToutt,3);
Elemcott(:,1)=-(vectort(:,3).*vectort(:,5)+vectort(:,4).*vectort(:,6))./(2*areat(:));
Elemcott(:,2)=-(vectort(:,1).*vectort(:,5)+vectort(:,2).*vectort(:,6))./(2*areat(:));
Elemcott(:,3)=-(vectort(:,3).*vectort(:,1)+vectort(:,4).*vectort(:,2))./(2*areat(:));



%% set boundary condition
%periodic for phi
[elem2dof,edge,~] = dofP2(elem);
elem2edge = elem2dof(:,4:6);
%[leftNode,rightNode,leftEdge,rightEdge] = getPerbd(node,edge);
%periodic for p,n and u
[outelem2dof,outedge,~] = dofP2(outelem);
outelem2edge = outelem2dof(:,4:6);
%[leftNodeS,rightNodeS,leftEdgeS,rightEdgeS] = getPerbd(outnode,outedge);


bdFlagPhi = setboundary(node, elem,'Dirichlet',' (abs(y-30)<0.00001)|(abs(abs(x)-30)<0.00001)');
%bdFlagPhi = setboundary(node, elem,'Dirichlet',' (abs(y-30)<0.00001)');
% figure(1)
% showmesh(node,elem);
% allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
% Dirichlet = allEdge((bdFlagPhi(:) == 1),:);
% Neumann = allEdge((bdFlagPhi(:) == 2),:);
% findedge(node,Dirichlet,[],'noindex','LineWidth',4,'Color','r');
% hold on;
% findedge(node,Neumann,[],'noindex','LineWidth',4,'Color','b');



bdFlagUz =  setboundary(outnode, outelem,'Dirichlet','(abs(abs(x)-30)<0.00001)|abs(y-30)<0.00001|(sqrt(x.^2+y.^2)-2)<0.00001');
bdFlagUr = setboundary(outnode, outelem,'Dirichlet','all');

% figure(2)
% showmesh(outnode,outelem);
% allEdge = [outelem(:,[2,3]); outelem(:,[3,1]); outelem(:,[1,2])];
% Dirichlet = allEdge((bdFlagUz(:) == 1),:);
% findedge(outnode,Dirichlet,[],'noindex','LineWidth',4,'Color','r');

bdFlagChar = setboundary(outnode,outelem,'Dirichlet','(abs(abs(x)-30)<0.00001)' );

figure(3)
showmesh(outnode,outelem);
allEdge = [outelem(:,[2,3]); outelem(:,[3,1]); outelem(:,[1,2])];
Dirichlet = allEdge((bdFlagChar(:) == 1),:);
findedge(outnode,Dirichlet,[],'noindex','LineWidth',4,'Color','r');





[outelem2doft,outedget,~] = dofP2(outelemt);
outelem2edget = outelem2doft(:,4:6);
%[leftNodeSt,rightNodeSt,leftEdgeSt,rightEdgeSt] = getPerbd(outnodet,outedget);
bdFlagt =  setboundary(outnodet,outelemt,'Neumann','abs(sqrt(x.^2+y.^2)-2.003)<0.00001');
bdFlagaux = setboundary(outnode,outelem,'Neumann','abs(sqrt(x.^2+y.^2)-2)<0.00001');







outallEdge = [outelem(:,[2,3]); outelem(:,[3,1]); outelem(:,[1,2])];
outNeumann = outallEdge((bdFlagaux(:) == 2),:);

outNeumann = sort(outNeumann,1);
% plot(outnode(outNeumann(:,1),1),outnode(outNeumann(:,1),2),'*');

% figure(4)
% showmesh(outnode,outelem);
% findedge(outnode,outNeumann,[],'noindex','LineWidth',4,'Color','r');

agn=pi-0.5*par.theta:-par.theta:0.5*par.theta;
sintheta=sin(agn');
costheta=cos(agn');

% T=auxstructure(outelemt);
% idxedge=T.elem2edge(:);
% idx=(bdFlagaux(:) == 2);
% edgeMark=idxedge(idx);
% elemMark=T.edge2elem(edgeMark,[1,3]);
%
elemMark=zeros(size(outNeumann,1),1);
for i=1:size(outNeumann,1)
    Lia=ismember(outelem,outNeumann(i,:));
    elemMark(i,1)=find(sum(Lia,2)==2);
    
end
showmesh(outnode,outelem(elemMark,:))

% for i=1:size(outNeumannt,1)
%     test(i,1)=ismember(outNeumannt(i,1),outelemt(elemMark(i,1),:));
%       test(i,2)=ismember(outNeumannt(i,2),outelemt(elemMark(i,1),:));
%
% end


%% pde
% get the Psi_c
pde = Pnpdata(SpotOld,par);

% showsolution(node,elem,pde.DPsic(node))





%% auxMesh
auxMesh = struct('node',node,'elem',elem,'outnode',outnode,'outelem',outelem,...
    'outnodeidx',outnodeidx,'outelemidx',outelemidx,'innernode',innernode,'innerelem',innerelem,...
    'innernodeidx',innernodeidx,'innerelemidx',innerelemidx,'elem2dof',elem2dof,'edge',edge,...
    'elem2edge',elem2edge,...
    'outelem2edge',outelem2edge,'outedge',outedge,'outelem2dof',outelem2dof,...
    'bdFlagPhi',bdFlagPhi,'bdFlagChar',bdFlagChar,'bdFlagUz',bdFlagUz,'bdFlagUr',bdFlagUr,'Elemcot',Elemcot);
auxMesht = struct('outnodet',outnodet,'outelemt',outelemt,...
    'outnodeidxt',outnodeidxt,'outelemidxt',outelemidxt,...
    'outelem2edget',outelem2edget,'outedget',outedget,'outelem2doft',outelem2doft,...
    'bdFlagt',bdFlagt,'Elemcott',Elemcott);
auxSphere = struct('elemMark',elemMark,'outNeumann',outNeumann, 'sintheta',sintheta,'costheta',costheta, 'theta',par.theta);
%% Assemble matrix
heateqMat=asmheatMat(auxMesh,auxMesht);
PoissMat=asmPoiMat(auxMesh,pde);
StokesMat = asmStokesP2P1(auxMesh,pde);
save auxMesh.mat auxMesh;
save auxMesht.mat auxMesht;
save pdedata.mat pde;
save auxSphere.mat auxSphere;
save par.mat par
save PoissMat.mat PoissMat;
save StokesMat.mat StokesMat;
save heateqMat.mat heateqMat;





