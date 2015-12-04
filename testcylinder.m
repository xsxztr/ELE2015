clear all
clc
load ('auxMesh.mat');
load ('pdedata.mat');
load('auxMesht.mat');
load('auxSphere.mat');
load('par.mat');
load('PoissMat.mat');
load('heateqMat.mat');
outnode = auxMesh.outnode;
outnodet = auxMesht.outnodet;
outelem = auxMesh.outelem;
outelemt = auxMesht.outelemt;
outnodeidx = auxMesh.outnodeidx;
Nout = size(outnode,1);
Nt= size(outnodet,1);
f1=zeros(Nout,1);%sin(outnode(:,1)/par.L*pi);
f2=zeros(Nt,1);%cos(outnodet(:,1)/par.L*pi);
node=auxMesh.node;
elem=auxMesh.elem;
surfchar=par.surchar;
%u=Cylinder(auxMesh,auxMesht,auxSphere,heateqMat,PoissMat,surfchar,f1,f2,pde);
phi=CylinderD(auxMesh,auxMesht,auxSphere,PoissMat,surfchar,f1,f2,pde);
u=zeros(Nout,2);
Iidx=1;
charge=getCharge(auxMesh,pde,phi,u,Iidx);

  StokesMat = asmStokesP2P1(auxMesh,pde);
 V1 = gradu(outnode,outelem,phi(outnodeidx));
 net2=zeros(Nout,1);
 V2=zeros(size(outelemt,1),2);
 
 svel=1;
 [u1,u2,p] = StokesP2P1(auxMesh,auxMesht,StokesMat,pde,charge,V1,net2,V2,svel);
  showsolution(outnode,outelem,u1);