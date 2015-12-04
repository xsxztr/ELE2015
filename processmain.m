clear all
clc
load ('auxMesh.mat');
load ('pdedata.mat');
load('auxMesht.mat');
load('auxSphere.mat');
load('par.mat');
load('PoissMat.mat');
load('heateqMat.mat');
load('StokesMat.mat');


outnode = auxMesh.outnode;
outnodet = auxMesht.outnodet;
node = auxMesh.node;
Nout = size(outnode,1);
Nt= size(outnodet,1);

Vmax = par.Vmax;
L = par.L;
den=1;
surfchar=1;
svel=1;

posOld = den*ones(Nout,1);
negOld = den*ones(Nout,1);
potentOld =  -Vmax*node(:,1)/(2*L);
option.nlItertol=0.000001;





[pos,neg,potent,u1,u2,press] = nlIter(potentOld,posOld,negOld,auxMesh,auxMesht,auxSphere,...
    pde,PoissMat,StokesMat,surfchar,svel,option);


fid1=fopen(['Vsur',sprintf('%i',sur),'n',sprintf('%i',den),'.txt'],'w');
fprintf(fid1,'%12.12i\r\n',potent);
fid2=fopen(['psur',sprintf('%i',sur),'n',sprintf('%i',den),'.txt'],'w');
fprintf(fid2,'%12.12i\r\n',pos);
fid3=fopen(['nsur',sprintf('%i',sur),'n',sprintf('%i',den),'.txt'],'w');
fprintf(fid3,'%12.12i\r\n',neg);
fid4=fopen(['U1sur',sprintf('%i',sur),'n',sprintf('%i',den),'.txt'],'w');
fprintf(fid4,'%12.12i\r\n',u1);
fid5=fopen(['U2sur',sprintf('%i',sur),'n',sprintf('%i',den),'.txt'],'w');
fprintf(fid5,'%12.12i\r\n',u2);
fid6=fopen(['presssur',sprintf('%i',sur),'n',sprintf('%i',den),'.txt'],'w');
fprintf(fid6,'%12.12i\r\n',press);


fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); fclose(fid5);fclose(fid6);
 