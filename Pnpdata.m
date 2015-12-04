function pdePnp = Pnpdata(Spot,par)



dy=par.dy; R=par.R;
R0=par.R0; L=par.L;

phi0=par.phi0;

C1=par.C1; C2=par.C2;

Dn=par.Dn; Dp=par.Dp;
Dk=par.Dk; Dcl=par.Dcl;

eta=par.eta;

epss=par.epss; epsf=par.epsf;
nH=par.nH; nOH=par.nOH;
nK=par.nK; ncl=par.ncl;

Vmax=par.Vmax;





%% ========================================================================
%% Homogenious Direchlet boundary condition for potent on r=R
    function z = phig_D(p)
        x = p(:,1);
        %z =zeros(size(x,1),1);
        z = -Vmax*x/(2*L);
    end
%% Homogenious Direchlet boundary condition for potent on r=R
    function z = phig_N(p)
        x = p(:,1);
        y = p(:,2);
        z =ones(size(x,1),1);
%         idx=(abs(y)<0.000001);
%         z(idx)=0;
        %z = -Vmax*x/(2*L);
    end
%%
   function z = charg_D(p)
        x = p(:,1);
        %z =zeros(size(x,1),1);
        z = ones(size(x,1),1);
   end

%% Homegenious direclet boundary condition for velocity ur on r=R and r=0
    function z = Urg_D(p)
        x=p(:,1);
        z=zeros(length(x),1);
    end
%% Homegenious direclet boundary condition for velocity uz on r=R
    function z = Uzg_D(p)
        x=p(:,1);
        y=p(:,2);
        z=ones(length(x),1);
        idx=(abs(sqrt(x.^2+y.^2)-2)<0.000001);
        z(idx)=0;
    end

%% trap function
    function z = Psic(p)
        x = p(:,1); y = p(:,2);
        r=sqrt((x-Spot).^2+y.^2);
        idx=(r-R<=dy);
        z=zeros(size(x,1),1);
        z(idx)=(1+cos(((r(idx)-2)*pi)/dy))/2*phi0;
    end
%% the gradient of the trap function
    function z = DPsic(p)
        x = p(:,1); y = p(:,2);
        r=sqrt((x-Spot).^2+y.^2);
        idx=(r-R<=dy);
        z=zeros(size(x,1),2);
        z(idx,1)=-sin(pi*(r(idx)-2)/dy)*pi/(2*dy)*phi0.*x(idx)./r(idx);
        z(idx,2)=-sin(pi*(r(idx)-2)/dy)*pi/(2*dy)*phi0.*y(idx)./r(idx);
    end



pdePnp = struct('c',1/10,'eta',eta,'C1',C1, 'C2',C2,'Dn',Dn,'Dp',Dp,...
    'Urg_D',@Urg_D,'Uzg_D',@Uzg_D,'phig_D',@phig_D,'phig_N',@phig_N,...
    'charg_D',@charg_D,'DPsic', @DPsic,'Psic',@Psic,'epss',epss,'epsf',epsf,...
    'nH',nH,'nOH',nOH,'nK',nK,'ncl',ncl,'Vmax',Vmax,'R',R,...
    'R0',R0,'dy',dy,'Dk',Dk,'Dcl',Dcl);
end