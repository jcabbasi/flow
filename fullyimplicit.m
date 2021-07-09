% In the name of God
% 2Phase 2dimension Water & Oil problem
% no flow Boundary
% Version : fuly implicit
clear
clc
% reservoir geometry and initial information
[XLength,YLength,ZLength,Phi,muw,muo,Cr,Co,Cw,initP,injQ,bhp,swi,sor,Boi,Bwi,betta_c,alpha_c,rw]=SystemInformation();

% some initial calculations
[delx,dely,delz,delt,Ax,Ay,Vb,Nx,Ny,Nt,InitP,bhp,Phi,XLength,YLength,ProdWells_X,ProdWells_Y,Inj_Wells_X,...
    Inj_Wells_Y,Perm]=initial_calculation(XLength,YLength,ZLength,initP,bhp,Phi);
dt=delt;

% Preallocations
POil=zeros(Ny,Nx,Nt); %Oil phase pressure
PWater=zeros(Ny,Nx,Nt); % Water phase pressure
sw=zeros(Ny,Nx,Nt); %water saturation
sw(:,:,1)=swi;
POil(:,:,1)=InitP; % insert inital pressure in Oil phase pressure in time step 1
[Pcow(:,:,1)]=PC_Func(sw(:,:,1),swi,sor); % calculate capillary pressure for time step 1
PWater(:,:,1)=InitP-Pcow(:,:,1); %  calculate Water phase pressyre in time step 1
coefficients=zeros((Ny)*(Nx)); % preallocation for coefficient matrix
constants=zeros(Ny*Nx,1); % preallocation for constants vector
bhpinj=zeros(Nt,1); % preallocation for bottom hole pressure for injection well
bhpinj(1,1)=3600;
Q_O=zeros(2,Nt); % preallocation for oil rate
Q_W=zeros(2,Nt); % preallocation for Water rate
f_w=zeros(2,Nt); % preallocation for water cut
X=zeros(Nx*Ny,1);
k=0;
for i=1:Ny
    for j=1:Nx
        X(k+1)=sw(i,j,1);
        X(k+2)=POil(i,j,1);
        k=k+2;
    end
end
n=2;
[G]=geometryfun(Nx,Ny,delx,dely,betta_c,Perm,Ax,Ay);

SW(:,:,1)=sw(:,:,n-1);
PW(:,:,1)=PWater(:,:,n-1);
PO(:,:,1)=POil(:,:,n-1);
%
dpp=1.001;dpm=0.999;dxp=0.002;
dsp=1.001;dsm=0.999;dxs=0.002;
%
err=10000;
v=1;
%%
injqw=zeros(Ny,Nx);
injqw(Inj_Wells_Y,Inj_Wells_X)=injQ;
injqo=zeros(Ny,Nx);
%
req=sqrt(delx.*dely/pi); %equivalent radius is calculated by peceman model
wc=2*pi*betta_c.*Perm.*delz./log(req/rw);  %productivity index
WC=zeros(Ny,Nx); % set productivity index for all grids equal zero
x=length(ProdWells_X);

for i=1:x
    WC(ProdWells_Y(i),ProdWells_X(i))=wc(ProdWells_Y(i),ProdWells_X(i));  % set productivity index for  grids which contains  production well
end
WC_inj(Inj_Wells_Y,Inj_Wells_X)=wc(Inj_Wells_Y,Inj_Wells_X);
n=2;
maxv=70;
maxdt=50;
dtchange=0.7;
tic
while n<10000
    SW(:,:,1)=sw(:,:,n-1);
    PW(:,:,1)=PWater(:,:,n-1);
    PO(:,:,1)=POil(:,:,n-1);
    err=zeros;
    err(1)=10000;
    v=1;
    itr=1;
    while err>10^-4
        if itr>maxv
            delt=delt*dtchange;
            itr=0;
        end
%         if v>maxv && err(v)<0.001
%             break;
%         end
        v=v+1;
        itr=itr+1;
        %
        Pv(:,:,1)=PO(:,:,v-1);
        Pv(:,:,2)=dpp*Pv(:,:,1);
        Pv(:,:,3)=dpm*Pv(:,:,1);
        %
        Sv(:,:,1)=SW(:,:,v-1);
        Sv(:,:,2)=dsp*Sv(:,:,1);
        Sv(:,:,3)=dsm*Sv(:,:,1);
        %
        [Pcow]=PC_Func(Sv(:,:,1),swi,sor);
        
        %
        Pwv(:,:,1)=PO(:,:,v-1)-Pcow;
        Pwv(:,:,2)=dpp*Pwv(:,:,1);
        Pwv(:,:,3)=dpm*Pwv(:,:,1);
        [Bo(:,:,1),Bw(:,:,1)]=FvfFunction(Pv(:,:,1),Pwv(:,:,1),Co,Cw,Ny,Nx);
        
        %
        [Kro(:,:,1),Krw(:,:,1)]=KR_Func(Sv(:,:,1),swi,sor);
        %
        for i=1:Ny
            for j=1:Nx
                GN=j+(i-1)*Nx;
                [fpo(GN,:,1)]=pressurefunterm(i,j,Nx,Ny,Pv(:,:,1),muo,Bo(:,:,1));
                %
                [fpw(GN,:,1)]=pressurefunterm(i,j,Nx,Ny,Pv(:,:,1),muw,Bw(:,:,1));
                %
                [fso(GN,:,1)]=saturationfun(i,j,Nx,Ny,Pv(:,:,1),Kro(:,:,1));
                %
                [fsw(GN,:,1)]=saturationfun(i,j,Nx,Ny,Pv(:,:,1),Krw(:,:,1));
                
            end
        end
        %
        [Cop(:,:,1),Cow(:,:,1),Cwp(:,:,1),Cww(:,:,1)]=TimeCoefficient(Vb,Phi,Sv(:,:,1),Boi,Bwi,Bo(:,:,1),Bw(:,:,1),Co,Cw,alpha_c,delt);
        
        %
        
        %%
        for i=1:Ny
            for j=1:Nx
                GN=j+(i-1)*Nx;
                To=G.*fpo(:,:,1).*fso(:,:,1);
                Tw=G.*fpw(:,:,1).*fsw(:,:,1);
                landaw=Krw(i,j)./muw./Bw(i,j);
                landao=Kro(i,j)./muo./Bo(i,j);
                Rvw=residualterm(i,j,Nx,Ny,Tw,Pv(:,:,1),Cwp(:,:,1),POil(:,:,n-1),Cww(:,:,1),Sv(:,:,1),sw(:,:,n-1),injqw,WC,landaw,bhp);
                Rvo=residualterm(i,j,Nx,Ny,To,Pv(:,:,1),Cop(:,:,1),POil(:,:,n-1),Cow(:,:,1),Sv(:,:,1),sw(:,:,n-1),injqo,WC,landao,bhp);
                %
                if i==1 && j==1
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    % i ,j+1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j+1)=Sv(i,j+1,2);
                    Svb(i,j+1)=Sv(i,j+1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j+1)=Pv(i,j+1,2);
                    Pvb(i,j+1)=Pv(i,j+1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j+1)=Pwv(i,j+1,2);
                    Pwvb(i,j+1)=Pwv(i,j+1,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j+1),Bwf(i,j+1)]=FvfFunction(Pvf(i,j+1),Pwvf(i,j+1),Co,Cw,Ny,Nx);
                    [Bob(i,j+1),Bwb(i,j+1)]=FvfFunction(Pvb(i,j+1),Pwvb(i,j+1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j+1),Krwf(i,j+1)]=KR_Func(Svf(i,j+1),swi,sor);
                    [Krob(i,j+1),Krwb(i,j+1)]=KR_Func(Svb(i,j+1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+1)-1:2*(GN+1))=[(Rwsf-Rwsb)/dxs/Sv(i,j+1,1) (Rwpf-Rwpb)/dxp/Pv(i,j+1,1);(Rosf-Rosb)/dxs/Sv(i,j+1,1) (Ropf-Ropb)/dxp/Pv(i,j+1,1)];
                    
                    % i+1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i+1,j)=Sv(i+1,j,2);
                    Svb(i+1,j)=Sv(i+1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i+1,j)=Pv(i+1,j,2);
                    Pvb(i+1,j)=Pv(i+1,j,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i+1,j)=Pwv(i+1,j,2);
                    Pwvb(i+1,j)=Pwv(i+1,j,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i+1,j),Bwf(i+1,j)]=FvfFunction(Pvf(i+1,j),Pwvf(i+1,j),Co,Cw,Ny,Nx);
                    [Bob(i+1,j),Bwb(i+1,j)]=FvfFunction(Pvb(i+1,j),Pwvb(i+1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i+1,j),Krwf(i+1,j)]=KR_Func(Svf(i+1,j),swi,sor);
                    [Krob(i+1,j),Krwb(i+1,j)]=KR_Func(Svb(i+1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+Nx)-1:2*(GN+Nx))=[(Rwsf-Rwsb)/dxs/Sv(i+1,j,1) (Rwpf-Rwpb)/dxp/Pv(i+1,j,1);(Rosf-Rosb)/dxs/Sv(i+1,j,1) (Ropf-Ropb)/dxp/Pv(i+1,j,1)];
                elseif i==1 && j==Nx
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    % i,j-1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j-1)=Sv(i,j-1,2);
                    Svb(i,j-1)=Sv(i,j-1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j-1)=Pv(i,j-1,2);
                    Pvb(i,j-1)=Pv(i,j-1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j-1)=Pwv(i,j-1,2);
                    Pwvb(i,j-1)=Pwv(i,j-1,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j-1),Bwf(i,j-1)]=FvfFunction(Pvf(i,j-1),Pwvf(i,j-1),Co,Cw,Ny,Nx);
                    [Bob(i,j-1),Bwb(i,j-1)]=FvfFunction(Pvb(i,j-1),Pwvb(i,j-1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j-1),Krwf(i,j-1)]=KR_Func(Svf(i,j-1),swi,sor);
                    [Krob(i,j-1),Krwb(i,j-1)]=KR_Func(Svb(i,j-1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN-1)-1:2*(GN-1))=[(Rwsf-Rwsb)/dxs/Sv(i,j-1,1) (Rwpf-Rwpb)/dxp/Pv(i,j-1,1);(Rosf-Rosb)/dxs/Sv(i,j-1,1) (Ropf-Ropb)/dxp/Pv(i,j-1,1)];
                    
                    % i+1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i+1,j)=Sv(i+1,j,2);
                    Svb(i+1,j)=Sv(i+1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i+1,j)=Pv(i+1,j,2);
                    Pvb(i+1,j)=Pv(i+1,j,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i+1,j)=Pwv(i+1,j,2);
                    Pwvb(i+1,j)=Pwv(i+1,j,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i+1,j),Bwf(i+1,j)]=FvfFunction(Pvf(i+1,j),Pwvf(i+1,j),Co,Cw,Ny,Nx);
                    [Bob(i+1,j),Bwb(i+1,j)]=FvfFunction(Pvb(i+1,j),Pwvb(i+1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i+1,j),Krwf(i+1,j)]=KR_Func(Svf(i+1,j),swi,sor);
                    [Krob(i+1,j),Krwb(i+1,j)]=KR_Func(Svb(i+1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+Nx)-1:2*(GN+Nx))=[(Rwsf-Rwsb)/dxs/Sv(i+1,j,1) (Rwpf-Rwpb)/dxp/Pv(i+1,j,1);(Rosf-Rosb)/dxs/Sv(i+1,j,1) (Ropf-Ropb)/dxp/Pv(i+1,j,1)];
                elseif i==Ny && j==1
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i-1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i-1,j)=Sv(i-1,j,2);
                    Svb(i-1,j)=Sv(i-1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i-1,j)=Pv(i-1,j,2);
                    Pvb(i-1,j)=Pv(i-1,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i-1,j)=Pwv(i-1,j,2);
                    Pwvb(i-1,j)=Pwv(i-1,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i-1,j),Bwf(i-1,j)]=FvfFunction(Pvf(i-1,j),Pwvf(i-1,j),Co,Cw,Ny,Nx);
                    [Bob(i-1,j),Bwb(i-1,j)]=FvfFunction(Pvb(i-1,j),Pwvb(i-1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i-1,j),Krwf(i-1,j)]=KR_Func(Svf(i-1,j),swi,sor);
                    [Krob(i-1,j),Krwb(i-1,j)]=KR_Func(Svb(i-1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*(GN-Nx)-1:2*(GN-Nx))=[(Rwsf-Rwsb)/dxs/Sv(i-1,j,1) (Rwpf-Rwpb)/dxp/Pv(i-1,j,1);(Rosf-Rosb)/dxs/Sv(i-1,j,1) (Ropf-Ropb)/dxp/Pv(i-1,j,1)];
                    % i ,j+1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j+1)=Sv(i,j+1,2);
                    Svb(i,j+1)=Sv(i,j+1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j+1)=Pv(i,j+1,2);
                    Pvb(i,j+1)=Pv(i,j+1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j+1)=Pwv(i,j+1,2);
                    Pwvb(i,j+1)=Pwv(i,j+1,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j+1),Bwf(i,j+1)]=FvfFunction(Pvf(i,j+1),Pwvf(i,j+1),Co,Cw,Ny,Nx);
                    [Bob(i,j+1),Bwb(i,j+1)]=FvfFunction(Pvb(i,j+1),Pwvb(i,j+1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j+1),Krwf(i,j+1)]=KR_Func(Svf(i,j+1),swi,sor);
                    [Krob(i,j+1),Krwb(i,j+1)]=KR_Func(Svb(i,j+1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+1)-1:2*(GN+1))=[(Rwsf-Rwsb)/dxs/Sv(i,j+1,1) (Rwpf-Rwpb)/dxp/Pv(i,j+1,1);(Rosf-Rosb)/dxs/Sv(i,j+1,1) (Ropf-Ropb)/dxp/Pv(i,j+1,1)];
                    
                elseif i==Ny && j==Nx
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i-1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i-1,j)=Sv(i-1,j,2);
                    Svb(i-1,j)=Sv(i-1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i-1,j)=Pv(i-1,j,2);
                    Pvb(i-1,j)=Pv(i-1,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i-1,j)=Pwv(i-1,j,2);
                    Pwvb(i-1,j)=Pwv(i-1,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i-1,j),Bwf(i-1,j)]=FvfFunction(Pvf(i-1,j),Pwvf(i-1,j),Co,Cw,Ny,Nx);
                    [Bob(i-1,j),Bwb(i-1,j)]=FvfFunction(Pvb(i-1,j),Pwvb(i-1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i-1,j),Krwf(i-1,j)]=KR_Func(Svf(i-1,j),swi,sor);
                    [Krob(i-1,j),Krwb(i-1,j)]=KR_Func(Svb(i-1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*(GN-Nx)-1:2*(GN-Nx))=[(Rwsf-Rwsb)/dxs/Sv(i-1,j,1) (Rwpf-Rwpb)/dxp/Pv(i-1,j,1);(Rosf-Rosb)/dxs/Sv(i-1,j,1) (Ropf-Ropb)/dxp/Pv(i-1,j,1)];
                    
                    % i,j-1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j-1)=Sv(i,j-1,2);
                    Svb(i,j-1)=Sv(i,j-1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j-1)=Pv(i,j-1,2);
                    Pvb(i,j-1)=Pv(i,j-1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j-1)=Pwv(i,j-1,2);
                    Pwvb(i,j-1)=Pwv(i,j-1,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j-1),Bwf(i,j-1)]=FvfFunction(Pvf(i,j-1),Pwvf(i,j-1),Co,Cw,Ny,Nx);
                    [Bob(i,j-1),Bwb(i,j-1)]=FvfFunction(Pvb(i,j-1),Pwvb(i,j-1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j-1),Krwf(i,j-1)]=KR_Func(Svf(i,j-1),swi,sor);
                    [Krob(i,j-1),Krwb(i,j-1)]=KR_Func(Svb(i,j-1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN-1)-1:2*(GN-1))=[(Rwsf-Rwsb)/dxs/Sv(i,j-1,1) (Rwpf-Rwpb)/dxp/Pv(i,j-1,1);(Rosf-Rosb)/dxs/Sv(i,j-1,1) (Ropf-Ropb)/dxp/Pv(i,j-1,1)];
                elseif i==1
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i,j-1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j-1)=Sv(i,j-1,2);
                    Svb(i,j-1)=Sv(i,j-1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j-1)=Pv(i,j-1,2);
                    Pvb(i,j-1)=Pv(i,j-1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j-1)=Pwv(i,j-1,2);
                    Pwvb(i,j-1)=Pwv(i,j-1,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j-1),Bwf(i,j-1)]=FvfFunction(Pvf(i,j-1),Pwvf(i,j-1),Co,Cw,Ny,Nx);
                    [Bob(i,j-1),Bwb(i,j-1)]=FvfFunction(Pvb(i,j-1),Pwvb(i,j-1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j-1),Krwf(i,j-1)]=KR_Func(Svf(i,j-1),swi,sor);
                    [Krob(i,j-1),Krwb(i,j-1)]=KR_Func(Svb(i,j-1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN-1)-1:2*(GN-1))=[(Rwsf-Rwsb)/dxs/Sv(i,j-1,1) (Rwpf-Rwpb)/dxp/Pv(i,j-1,1);(Rosf-Rosb)/dxs/Sv(i,j-1,1) (Ropf-Ropb)/dxp/Pv(i,j-1,1)];
                    
                    % i ,j+1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j+1)=Sv(i,j+1,2);
                    Svb(i,j+1)=Sv(i,j+1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j+1)=Pv(i,j+1,2);
                    Pvb(i,j+1)=Pv(i,j+1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j+1)=Pwv(i,j+1,2);
                    Pwvb(i,j+1)=Pwv(i,j+1,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j+1),Bwf(i,j+1)]=FvfFunction(Pvf(i,j+1),Pwvf(i,j+1),Co,Cw,Ny,Nx);
                    [Bob(i,j+1),Bwb(i,j+1)]=FvfFunction(Pvb(i,j+1),Pwvb(i,j+1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j+1),Krwf(i,j+1)]=KR_Func(Svf(i,j+1),swi,sor);
                    [Krob(i,j+1),Krwb(i,j+1)]=KR_Func(Svb(i,j+1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+1)-1:2*(GN+1))=[(Rwsf-Rwsb)/dxs/Sv(i,j+1,1) (Rwpf-Rwpb)/dxp/Pv(i,j+1,1);(Rosf-Rosb)/dxs/Sv(i,j+1,1) (Ropf-Ropb)/dxp/Pv(i,j+1,1)];
                    
                    % i+1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i+1,j)=Sv(i+1,j,2);
                    Svb(i+1,j)=Sv(i+1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i+1,j)=Pv(i+1,j,2);
                    Pvb(i+1,j)=Pv(i+1,j,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i+1,j)=Pwv(i+1,j,2);
                    Pwvb(i+1,j)=Pwv(i+1,j,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i+1,j),Bwf(i+1,j)]=FvfFunction(Pvf(i+1,j),Pwvf(i+1,j),Co,Cw,Ny,Nx);
                    [Bob(i+1,j),Bwb(i+1,j)]=FvfFunction(Pvb(i+1,j),Pwvb(i+1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i+1,j),Krwf(i+1,j)]=KR_Func(Svf(i+1,j),swi,sor);
                    [Krob(i+1,j),Krwb(i+1,j)]=KR_Func(Svb(i+1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+Nx)-1:2*(GN+Nx))=[(Rwsf-Rwsb)/dxs/Sv(i+1,j,1) (Rwpf-Rwpb)/dxp/Pv(i+1,j,1);(Rosf-Rosb)/dxs/Sv(i+1,j,1) (Ropf-Ropb)/dxp/Pv(i+1,j,1)];
                elseif j==1
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i-1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i-1,j)=Sv(i-1,j,2);
                    Svb(i-1,j)=Sv(i-1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i-1,j)=Pv(i-1,j,2);
                    Pvb(i-1,j)=Pv(i-1,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i-1,j)=Pwv(i-1,j,2);
                    Pwvb(i-1,j)=Pwv(i-1,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i-1,j),Bwf(i-1,j)]=FvfFunction(Pvf(i-1,j),Pwvf(i-1,j),Co,Cw,Ny,Nx);
                    [Bob(i-1,j),Bwb(i-1,j)]=FvfFunction(Pvb(i-1,j),Pwvb(i-1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i-1,j),Krwf(i-1,j)]=KR_Func(Svf(i-1,j),swi,sor);
                    [Krob(i-1,j),Krwb(i-1,j)]=KR_Func(Svb(i-1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*(GN-Nx)-1:2*(GN-Nx))=[(Rwsf-Rwsb)/dxs/Sv(i-1,j,1) (Rwpf-Rwpb)/dxp/Pv(i-1,j,1);(Rosf-Rosb)/dxs/Sv(i-1,j,1) (Ropf-Ropb)/dxp/Pv(i-1,j,1)];
                    
                    % i ,j+1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j+1)=Sv(i,j+1,2);
                    Svb(i,j+1)=Sv(i,j+1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j+1)=Pv(i,j+1,2);
                    Pvb(i,j+1)=Pv(i,j+1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j+1)=Pwv(i,j+1,2);
                    Pwvb(i,j+1)=Pwv(i,j+1,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j+1),Bwf(i,j+1)]=FvfFunction(Pvf(i,j+1),Pwvf(i,j+1),Co,Cw,Ny,Nx);
                    [Bob(i,j+1),Bwb(i,j+1)]=FvfFunction(Pvb(i,j+1),Pwvb(i,j+1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j+1),Krwf(i,j+1)]=KR_Func(Svf(i,j+1),swi,sor);
                    [Krob(i,j+1),Krwb(i,j+1)]=KR_Func(Svb(i,j+1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+1)-1:2*(GN+1))=[(Rwsf-Rwsb)/dxs/Sv(i,j+1,1) (Rwpf-Rwpb)/dxp/Pv(i,j+1,1);(Rosf-Rosb)/dxs/Sv(i,j+1,1) (Ropf-Ropb)/dxp/Pv(i,j+1,1)];
                    
                    % i+1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i+1,j)=Sv(i+1,j,2);
                    Svb(i+1,j)=Sv(i+1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i+1,j)=Pv(i+1,j,2);
                    Pvb(i+1,j)=Pv(i+1,j,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i+1,j)=Pwv(i+1,j,2);
                    Pwvb(i+1,j)=Pwv(i+1,j,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i+1,j),Bwf(i+1,j)]=FvfFunction(Pvf(i+1,j),Pwvf(i+1,j),Co,Cw,Ny,Nx);
                    [Bob(i+1,j),Bwb(i+1,j)]=FvfFunction(Pvb(i+1,j),Pwvb(i+1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i+1,j),Krwf(i+1,j)]=KR_Func(Svf(i+1,j),swi,sor);
                    [Krob(i+1,j),Krwb(i+1,j)]=KR_Func(Svb(i+1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+Nx)-1:2*(GN+Nx))=[(Rwsf-Rwsb)/dxs/Sv(i+1,j,1) (Rwpf-Rwpb)/dxp/Pv(i+1,j,1);(Rosf-Rosb)/dxs/Sv(i+1,j,1) (Ropf-Ropb)/dxp/Pv(i+1,j,1)];
                elseif i==Ny
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i-1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i-1,j)=Sv(i-1,j,2);
                    Svb(i-1,j)=Sv(i-1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i-1,j)=Pv(i-1,j,2);
                    Pvb(i-1,j)=Pv(i-1,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i-1,j)=Pwv(i-1,j,2);
                    Pwvb(i-1,j)=Pwv(i-1,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i-1,j),Bwf(i-1,j)]=FvfFunction(Pvf(i-1,j),Pwvf(i-1,j),Co,Cw,Ny,Nx);
                    [Bob(i-1,j),Bwb(i-1,j)]=FvfFunction(Pvb(i-1,j),Pwvb(i-1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i-1,j),Krwf(i-1,j)]=KR_Func(Svf(i-1,j),swi,sor);
                    [Krob(i-1,j),Krwb(i-1,j)]=KR_Func(Svb(i-1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*(GN-Nx)-1:2*(GN-Nx))=[(Rwsf-Rwsb)/dxs/Sv(i-1,j,1) (Rwpf-Rwpb)/dxp/Pv(i-1,j,1);(Rosf-Rosb)/dxs/Sv(i-1,j,1) (Ropf-Ropb)/dxp/Pv(i-1,j,1)];
                    
                    % i,j-1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j-1)=Sv(i,j-1,2);
                    Svb(i,j-1)=Sv(i,j-1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j-1)=Pv(i,j-1,2);
                    Pvb(i,j-1)=Pv(i,j-1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j-1)=Pwv(i,j-1,2);
                    Pwvb(i,j-1)=Pwv(i,j-1,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j-1),Bwf(i,j-1)]=FvfFunction(Pvf(i,j-1),Pwvf(i,j-1),Co,Cw,Ny,Nx);
                    [Bob(i,j-1),Bwb(i,j-1)]=FvfFunction(Pvb(i,j-1),Pwvb(i,j-1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j-1),Krwf(i,j-1)]=KR_Func(Svf(i,j-1),swi,sor);
                    [Krob(i,j-1),Krwb(i,j-1)]=KR_Func(Svb(i,j-1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN-1)-1:2*(GN-1))=[(Rwsf-Rwsb)/dxs/Sv(i,j-1,1) (Rwpf-Rwpb)/dxp/Pv(i,j-1,1);(Rosf-Rosb)/dxs/Sv(i,j-1,1) (Ropf-Ropb)/dxp/Pv(i,j-1,1)];
                    
                    % i ,j+1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j+1)=Sv(i,j+1,2);
                    Svb(i,j+1)=Sv(i,j+1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j+1)=Pv(i,j+1,2);
                    Pvb(i,j+1)=Pv(i,j+1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j+1)=Pwv(i,j+1,2);
                    Pwvb(i,j+1)=Pwv(i,j+1,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j+1),Bwf(i,j+1)]=FvfFunction(Pvf(i,j+1),Pwvf(i,j+1),Co,Cw,Ny,Nx);
                    [Bob(i,j+1),Bwb(i,j+1)]=FvfFunction(Pvb(i,j+1),Pwvb(i,j+1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j+1),Krwf(i,j+1)]=KR_Func(Svf(i,j+1),swi,sor);
                    [Krob(i,j+1),Krwb(i,j+1)]=KR_Func(Svb(i,j+1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+1)-1:2*(GN+1))=[(Rwsf-Rwsb)/dxs/Sv(i,j+1,1) (Rwpf-Rwpb)/dxp/Pv(i,j+1,1);(Rosf-Rosb)/dxs/Sv(i,j+1,1) (Ropf-Ropb)/dxp/Pv(i,j+1,1)];
                elseif j==Nx
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i-1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i-1,j)=Sv(i-1,j,2);
                    Svb(i-1,j)=Sv(i-1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i-1,j)=Pv(i-1,j,2);
                    Pvb(i-1,j)=Pv(i-1,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i-1,j)=Pwv(i-1,j,2);
                    Pwvb(i-1,j)=Pwv(i-1,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i-1,j),Bwf(i-1,j)]=FvfFunction(Pvf(i-1,j),Pwvf(i-1,j),Co,Cw,Ny,Nx);
                    [Bob(i-1,j),Bwb(i-1,j)]=FvfFunction(Pvb(i-1,j),Pwvb(i-1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i-1,j),Krwf(i-1,j)]=KR_Func(Svf(i-1,j),swi,sor);
                    [Krob(i-1,j),Krwb(i-1,j)]=KR_Func(Svb(i-1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*(GN-Nx)-1:2*(GN-Nx))=[(Rwsf-Rwsb)/dxs/Sv(i-1,j,1) (Rwpf-Rwpb)/dxp/Pv(i-1,j,1);(Rosf-Rosb)/dxs/Sv(i-1,j,1) (Ropf-Ropb)/dxp/Pv(i-1,j,1)];
                    
                    % i,j-1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j-1)=Sv(i,j-1,2);
                    Svb(i,j-1)=Sv(i,j-1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j-1)=Pv(i,j-1,2);
                    Pvb(i,j-1)=Pv(i,j-1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j-1)=Pwv(i,j-1,2);
                    Pwvb(i,j-1)=Pwv(i,j-1,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j-1),Bwf(i,j-1)]=FvfFunction(Pvf(i,j-1),Pwvf(i,j-1),Co,Cw,Ny,Nx);
                    [Bob(i,j-1),Bwb(i,j-1)]=FvfFunction(Pvb(i,j-1),Pwvb(i,j-1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j-1),Krwf(i,j-1)]=KR_Func(Svf(i,j-1),swi,sor);
                    [Krob(i,j-1),Krwb(i,j-1)]=KR_Func(Svb(i,j-1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN-1)-1:2*(GN-1))=[(Rwsf-Rwsb)/dxs/Sv(i,j-1,1) (Rwpf-Rwpb)/dxp/Pv(i,j-1,1);(Rosf-Rosb)/dxs/Sv(i,j-1,1) (Ropf-Ropb)/dxp/Pv(i,j-1,1)];
                    % i+1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i+1,j)=Sv(i+1,j,2);
                    Svb(i+1,j)=Sv(i+1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i+1,j)=Pv(i+1,j,2);
                    Pvb(i+1,j)=Pv(i+1,j,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i+1,j)=Pwv(i+1,j,2);
                    Pwvb(i+1,j)=Pwv(i+1,j,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i+1,j),Bwf(i+1,j)]=FvfFunction(Pvf(i+1,j),Pwvf(i+1,j),Co,Cw,Ny,Nx);
                    [Bob(i+1,j),Bwb(i+1,j)]=FvfFunction(Pvb(i+1,j),Pwvb(i+1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i+1,j),Krwf(i+1,j)]=KR_Func(Svf(i+1,j),swi,sor);
                    [Krob(i+1,j),Krwb(i+1,j)]=KR_Func(Svb(i+1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+Nx)-1:2*(GN+Nx))=[(Rwsf-Rwsb)/dxs/Sv(i+1,j,1) (Rwpf-Rwpb)/dxp/Pv(i+1,j,1);(Rosf-Rosb)/dxs/Sv(i+1,j,1) (Ropf-Ropb)/dxp/Pv(i+1,j,1)];
                else
                    % i,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j)=Sv(i,j,2);
                    Svb(i,j)=Sv(i,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j)=Pv(i,j,2);
                    Pvb(i,j)=Pv(i,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j)=Pwv(i,j,2);
                    Pwvb(i,j)=Pwv(i,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j),Bwf(i,j)]=FvfFunction(Pvf(i,j),Pwvf(i,j),Co,Cw,Ny,Nx);
                    [Bob(i,j),Bwb(i,j)]=FvfFunction(Pvb(i,j),Pwvb(i,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j),Krwf(i,j)]=KR_Func(Svf(i,j),swi,sor);
                    [Krob(i,j),Krwb(i,j)]=KR_Func(Svb(i,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bof(i,j);
                    landaopb=Kro(i,j)/muo/Bob(i,j);
                    landaosf=Krof(i,j)/muo/Bo(i,j);
                    landaosb=Krob(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bwf(i,j);
                    landawpb=Krw(i,j)/muw/Bwb(i,j);
                    landawsf=Krwf(i,j)/muw/Bw(i,j);
                    landawsb=Krwb(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*GN-1:2*GN)=[(Rwsf-Rwsb)/dxs/Sv(i,j,1) (Rwpf-Rwpb)/dxp/Pv(i,j,1);(Rosf-Rosb)/dxs/Sv(i,j,1) (Ropf-Ropb)/dxp/Pv(i,j,1)];
                    
                    % i-1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i-1,j)=Sv(i-1,j,2);
                    Svb(i-1,j)=Sv(i-1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i-1,j)=Pv(i-1,j,2);
                    Pvb(i-1,j)=Pv(i-1,j,3);
                    %
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i-1,j)=Pwv(i-1,j,2);
                    Pwvb(i-1,j)=Pwv(i-1,j,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i-1,j),Bwf(i-1,j)]=FvfFunction(Pvf(i-1,j),Pwvf(i-1,j),Co,Cw,Ny,Nx);
                    [Bob(i-1,j),Bwb(i-1,j)]=FvfFunction(Pvb(i-1,j),Pwvb(i-1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i-1,j),Krwf(i-1,j)]=KR_Func(Svf(i-1,j),swi,sor);
                    [Krob(i-1,j),Krwb(i-1,j)]=KR_Func(Svb(i-1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    
                    J(2*GN-1:2*GN,2*(GN-Nx)-1:2*(GN-Nx))=[(Rwsf-Rwsb)/dxs/Sv(i-1,j,1) (Rwpf-Rwpb)/dxp/Pv(i-1,j,1);(Rosf-Rosb)/dxs/Sv(i-1,j,1) (Ropf-Ropb)/dxp/Pv(i-1,j,1)];
                    
                    % i,j-1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j-1)=Sv(i,j-1,2);
                    Svb(i,j-1)=Sv(i,j-1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j-1)=Pv(i,j-1,2);
                    Pvb(i,j-1)=Pv(i,j-1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j-1)=Pwv(i,j-1,2);
                    Pwvb(i,j-1)=Pwv(i,j-1,3);
                    %
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j-1),Bwf(i,j-1)]=FvfFunction(Pvf(i,j-1),Pwvf(i,j-1),Co,Cw,Ny,Nx);
                    [Bob(i,j-1),Bwb(i,j-1)]=FvfFunction(Pvb(i,j-1),Pwvb(i,j-1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j-1),Krwf(i,j-1)]=KR_Func(Svf(i,j-1),swi,sor);
                    [Krob(i,j-1),Krwb(i,j-1)]=KR_Func(Svb(i,j-1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN-1)-1:2*(GN-1))=[(Rwsf-Rwsb)/dxs/Sv(i,j-1,1) (Rwpf-Rwpb)/dxp/Pv(i,j-1,1);(Rosf-Rosb)/dxs/Sv(i,j-1,1) (Ropf-Ropb)/dxp/Pv(i,j-1,1)];
                    
                    % i ,j+1
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i,j+1)=Sv(i,j+1,2);
                    Svb(i,j+1)=Sv(i,j+1,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i,j+1)=Pv(i,j+1,2);
                    Pvb(i,j+1)=Pv(i,j+1,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i,j+1)=Pwv(i,j+1,2);
                    Pwvb(i,j+1)=Pwv(i,j+1,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i,j+1),Bwf(i,j+1)]=FvfFunction(Pvf(i,j+1),Pwvf(i,j+1),Co,Cw,Ny,Nx);
                    [Bob(i,j+1),Bwb(i,j+1)]=FvfFunction(Pvb(i,j+1),Pwvb(i,j+1),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i,j+1),Krwf(i,j+1)]=KR_Func(Svf(i,j+1),swi,sor);
                    [Krob(i,j+1),Krwb(i,j+1)]=KR_Func(Svb(i,j+1),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+1)-1:2*(GN+1))=[(Rwsf-Rwsb)/dxs/Sv(i,j+1,1) (Rwpf-Rwpb)/dxp/Pv(i,j+1,1);(Rosf-Rosb)/dxs/Sv(i,j+1,1) (Ropf-Ropb)/dxp/Pv(i,j+1,1)];
                    
                    % i+1 ,j
                    Svf=Sv(:,:,1);Svb=Sv(:,:,1);
                    Svf(i+1,j)=Sv(i+1,j,2);
                    Svb(i+1,j)=Sv(i+1,j,3);
                    
                    Pvf=Pv(:,:,1);Pvb=Pv(:,:,1);
                    Pvf(i+1,j)=Pv(i+1,j,2);
                    Pvb(i+1,j)=Pv(i+1,j,3);
                    
                    Pwvf=Pwv(:,:,1);Pwvb=Pwv(:,:,1);
                    Pwvf(i+1,j)=Pwv(i+1,j,2);
                    Pwvb(i+1,j)=Pwv(i+1,j,3);
                    %
                    Bof=Bo(:,:,1);
                    Bob=Bo(:,:,1);
                    Bwf=Bw(:,:,1);
                    Bwb=Bw(:,:,1);
                    [Bof(i+1,j),Bwf(i+1,j)]=FvfFunction(Pvf(i+1,j),Pwvf(i+1,j),Co,Cw,Ny,Nx);
                    [Bob(i+1,j),Bwb(i+1,j)]=FvfFunction(Pvb(i+1,j),Pwvb(i+1,j),Co,Cw,Ny,Nx);
                    %
                    Krof=Kro(:,:,1);
                    Krob=Kro(:,:,1);
                    Krwf=Krw(:,:,1);
                    Krwb=Krw(:,:,1);
                    
                    [Krof(i+1,j),Krwf(i+1,j)]=KR_Func(Svf(i+1,j),swi,sor);
                    [Krob(i+1,j),Krwb(i+1,j)]=KR_Func(Svb(i+1,j),swi,sor);
                    %
                    [fppsof]=Trans_PS_Fun(Nx,Ny,Pvf,Kro,muo,Bof);
                    [fppsob]=Trans_PS_Fun(Nx,Ny,Pvb,Kro,muo,Bob);
                    [fppswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krw,muw,Bwf);
                    [fppswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krw,muw,Bwb);
                    %
                    Topf=G.*fppsof;
                    Topb=G.*fppsob;
                    
                    Twpf=G.*fppswf;
                    Twpb=G.*fppswb;
                    %
                    [fpssof]=Trans_PS_Fun(Nx,Ny,Pvf,Krof,muo,Bo);
                    [fpssob]=Trans_PS_Fun(Nx,Ny,Pvb,Krob,muo,Bo);
                    [fpsswf]=Trans_PS_Fun(Nx,Ny,Pvf,Krwf,muw,Bw);
                    [fpsswb]=Trans_PS_Fun(Nx,Ny,Pvb,Krwb,muw,Bw);
                    %
                    Tosf=G.*fpssof;
                    Tosb=G.*fpssob;
                    
                    Twsf=G.*fpsswf;
                    Twsb=G.*fpsswb;
                    %
                    [Copf,Cowf,Cwpf,Cwwf]=TimeCoefficient(Vb,Phi,Svf,Boi,Bwi,Bof,Bwf,Co,Cw,alpha_c,delt);
                    [Copb,Cowb,Cwpb,Cwwb]=TimeCoefficient(Vb,Phi,Svb,Boi,Bwi,Bob,Bwb,Co,Cw,alpha_c,delt);
                    %
                    landaopf=Kro(i,j)/muo/Bo(i,j);
                    landaopb=Kro(i,j)/muo/Bo(i,j);
                    landaosf=Kro(i,j)/muo/Bo(i,j);
                    landaosb=Kro(i,j)/muo/Bo(i,j);
                    
                    landawpf=Krw(i,j)/muw/Bw(i,j);
                    landawpb=Krw(i,j)/muw/Bw(i,j);
                    landawsf=Krw(i,j)/muw/Bw(i,j);
                    landawsb=Krw(i,j)/muw/Bw(i,j);
                    %
                    [Rwsf]=residualterm(i,j,Nx,Ny,Twsf,Pwv(:,:,1),Cwpf,POil(:,:,n-1),Cww(:,:,1),Svf,sw(:,:,n-1),injqw,WC,landawsf,bhp);
                    [Rwsb]=residualterm(i,j,Nx,Ny,Twsb,Pwv(:,:,1),Cwpb,POil(:,:,n-1),Cww(:,:,1),Svb,sw(:,:,n-1),injqw,WC,landawsb,bhp);
                    
                    [Rwpf]=residualterm(i,j,Nx,Ny,Twpf,Pwvf,Cwp(:,:,1),POil(:,:,n-1),Cwwf,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpf,bhp);
                    [Rwpb]=residualterm(i,j,Nx,Ny,Twpb,Pwvb,Cwp(:,:,1),POil(:,:,n-1),Cwwb,Sv(:,:,1),sw(:,:,n-1),injqw,WC,landawpb,bhp);
                    
                    [Rosf]=residualterm(i,j,Nx,Ny,Tosf,Pv(:,:,1),Copf,POil(:,:,n-1),Cow(:,:,1),Svf,sw(:,:,n-1),injqo,WC,landaosf,bhp);
                    [Rosb]=residualterm(i,j,Nx,Ny,Tosb,Pv(:,:,1),Copb,POil(:,:,n-1),Cow(:,:,1),Svb,sw(:,:,n-1),injqo,WC,landaosb,bhp);
                    
                    [Ropf]=residualterm(i,j,Nx,Ny,Topf,Pvf,Cop(:,:,1),POil(:,:,n-1),Cowf,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopf,bhp);
                    [Ropb]=residualterm(i,j,Nx,Ny,Topb,Pvb,Cop(:,:,1),POil(:,:,n-1),Cowb,Sv(:,:,1),sw(:,:,n-1),injqo,WC,landaopb,bhp);
                    J(2*GN-1:2*GN,2*(GN+Nx)-1:2*(GN+Nx))=[(Rwsf-Rwsb)/dxs/Sv(i+1,j,1) (Rwpf-Rwpb)/dxp/Pv(i+1,j,1);(Rosf-Rosb)/dxs/Sv(i+1,j,1) (Ropf-Ropb)/dxp/Pv(i+1,j,1)];
                    
                end
                
                Rv(2*GN-1:2*GN,1)=[Rvw;Rvo];
            end
        end
        
        DELX=J\-Rv;
        X=X+DELX;
        kk=0;
        for i=1:Ny
            for j=1:Nx
                SW(i,j,v)=X(kk+1);
                PO(i,j,v)=X(kk+2);
                PW(i,j,v)=X(kk+2);
                kk=kk+2;
            end
        end
        
        err(1,v)=sum(abs(DELX));
        err1(1,v)=sum(abs((Rv)));
        
    end
   % make time step vector
    DTT(1,n)=delt; 
    if (maxv-itr)>0
        delt=delt/dtchange;
        if delt<maxdt
        else
            delt=delt*dtchange;
        end
        
    end
    sw(:,:,n)=SW(:,:,v);
    PWater(:,:,n)=PW(:,:,v);
    POil(:,:,n)=PO(:,:,v);
    
    % Calculate injection well's bottom hole pressure.
    landao=mobilityn(Kro(Inj_Wells_Y,Inj_Wells_X),muo,Bo(Inj_Wells_Y,Inj_Wells_X));landaw=mobilityn(Krw(Inj_Wells_Y,Inj_Wells_X),muw,Bw(Inj_Wells_Y,Inj_Wells_X));
    bhpinj(n,1)=PWater(Inj_Wells_Y,Inj_Wells_X,n)+injQ/(WC_inj(Inj_Wells_Y,Inj_Wells_X)*((Bo(Inj_Wells_Y,Inj_Wells_X)/Bw(Inj_Wells_Y,Inj_Wells_X))*landao+landaw));
    %     % Flow rate calculation
    x=length(ProdWells_X);
    for i=1:x
        
        landao=mobilityn(Kro(ProdWells_Y(i),ProdWells_X(i)),muo,Bo(ProdWells_Y(i),ProdWells_X(i)));landaw=mobilityn(Krw(ProdWells_Y(i),ProdWells_X(i)),muw,Bw(ProdWells_Y(i),ProdWells_X(i)));
        
        Q_O(i,n)=WC(ProdWells_Y(i),ProdWells_X(i))*landao*(POil(ProdWells_Y(i),ProdWells_X(i),n)-bhp(ProdWells_Y(i),ProdWells_X(i)));
        Q_W(i,n)=WC(ProdWells_Y(i),ProdWells_X(i))*landaw*(PWater(ProdWells_Y(i),ProdWells_X(i),n)-bhp(ProdWells_Y(i),ProdWells_X(i)));
        
        %         % Water Cut Calculation.
        f_w(i,n)=Q_W(i,n)/(Q_W(i,n)+Q_O(i,n));
    end
    if f_w(1,n)>=0.95 && f_w(2,n) >=0.95 && f_w(3,n) >=0.95 && f_w(4,n) >=0.95  % condition for stop simulation
        break;
    end
    n=n+1;
end
toc
% Initial oil in place
OIP=sum(sum(Vb.*Phi.*(1-swi)./Boi))/alpha_c; % calcilate  initial oil in place

%  totally Water & Oil Daily & Cumulative Production.
WaterDaily_Production=zeros(1,n); %preallocation for daily  production
OilDaily_Production=zeros(1,n);
WaterCumulative_Production=zeros(1,n);%preallocation for cumulative  production
OilCumulative_Production=zeros(1,n);

% recovery factor
n=n-1;
R_F=zeros(1,n); % recovery facor preallocation
for i=2:n
    %  totally daily and cumulative productions 
    OilDaily_Production(1,i)=sum(Q_O(:,i)*DTT(1,i));
    OilCumulative_Production(1,i) =OilCumulative_Production(1,i-1)+sum(Q_O(:,i)*DTT(1,i));
    WaterDaily_Production(1,i)=sum(Q_W(:,i)*DTT(1,i));
    WaterCumulative_Production(1,i) =WaterCumulative_Production(1,i-1)+sum(Q_W(:,i)*DTT(1,i));
    % Recovery Factor
    R_F(1,i)= OilCumulative_Production(1,i)/OIP;
end

% preallocation production for well 1
Oildailyprod1=zeros(1,n);
Oilcumprod1=zeros(1,n);
%
waterdailyprod1=zeros(1,n);
watercumprod1=zeros(1,n);
% preallocation production for well 2
Oildailyprod2=zeros(1,n);
Oilcumprod2=zeros(1,n);
%
waterdailyprod2=zeros(1,n);
watercumprod2=zeros(1,n);
% preallocation production for well 3
Oildailyprod3=zeros(1,n);
Oilcumprod3=zeros(1,n);
%
waterdailyprod3=zeros(1,n);
watercumprod3=zeros(1,n);
% preallocation production for well 4
Oildailyprod4=zeros(1,n);
Oilcumprod4=zeros(1,n);
%
waterdailyprod4=zeros(1,n);
watercumprod4=zeros(1,n);
for i=2:n
    
    % daily and cumulative production for well num 1
    Oildailyprod1(1,i)=Q_O(1,i)*DTT(1,i);
    Oilcumprod1(1,i)=Oilcumprod1(1,i-1)+Q_O(1,i)*DTT(1,i);
    %
    waterdailyprod1(1,i)=Q_W(1,i)*DTT(1,i);
    watercumprod1(1,i)=watercumprod1(1,i-1)+Q_W(1,i)*DTT(1,i);
    % daily and cumulative production for well num 2
    Oildailyprod2(1,i)=Q_O(2,i)*DTT(1,i);
    Oilcumprod2(1,i)=Oilcumprod2(1,i-1)+Q_O(2,i)*DTT(1,i);
    %
    waterdailyprod2(1,i)=Q_W(2,i)*DTT(1,i);
    watercumprod2(1,i)=watercumprod2(1,i-1)+Q_W(2,i)*DTT(1,i);
    % daily and cumulative production for well num 3
    Oildailyprod3(1,i)=Q_O(3,i)*DTT(1,i);
    Oilcumprod3(1,i)=Oilcumprod3(1,i-1)+Q_O(3,i)*DTT(1,i);
    %
    waterdailyprod3(1,i)=Q_W(3,i)*DTT(1,i);
    watercumprod3(1,i)=watercumprod3(1,i-1)+Q_W(3,i)*DTT(1,i);
    % daily and cumulative production for well num 4
    Oildailyprod4(1,i)=Q_O(4,i)*DTT(1,i);
    Oilcumprod4(1,i)=Oilcumprod4(1,i-1)+Q_O(4,i)*DTT(1,i);
    %
    waterdailyprod4(1,i)=Q_W(4,i)*DTT(1,i);
    watercumprod4(1,i)=watercumprod4(1,i-1)+Q_W(4,i)*DTT(1,i);
    
end

% average Pressure and Water saturation in reservoir
AvPressure=zeros(1,n); %preallocation for Average Oil pressure
AvWaterSat=zeros(1,n);%preallocation for Average water pressure
for i=1:n
    AvPressure(1,i)=sum(sum(POil(:,:,i).*delx(:,:).*dely(:,:)))/(XLength*YLength);
    AvWaterSat(1,i)=sum(sum(sw(:,:,i).*delx(:,:).*dely(:,:)))/(XLength*YLength);
end

%%
% make time vector for figures(because of adaptive time step)
time=zeros;
for i=1:n
time(i,1)=sum(DTT(1,1:i));
end
%%
% Plot Injection well bottom hole pressure vs time
figure(1)
plot(time(2:n),bhpinj(2:n),'r.-')
hand=title(' Injection Well Bottom Hole Pressure');
set(hand,'fontsize',14);
xlabel('Time(days)','fontsize',13);
ylabel('pressure(psi)','fontsize',13);



% plot totally  daily and cumulative productions
% Water Production
figure(2)
title('water production for all wells')
yyaxis right
plot(time,WaterCumulative_Production(1:n),'r.-')
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,WaterDaily_Production(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);

% Oil Production
figure(3)
title('oil production for all wells')
yyaxis right
plot(time,OilCumulative_Production(1:n),'r.-')
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,OilDaily_Production(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%



% well 1: Oil Production 
figure(4)
title('oil production for well 1')
yyaxis right
plot(time,Oilcumprod1(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,Oildailyprod1(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% well 2: Oil Production 
figure(5)
title('oil production for well 2')
yyaxis right
plot(time,Oilcumprod2(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,Oildailyprod2(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% well 3: Oil Production 
figure(6)
title('oil production for well 3')
yyaxis right
plot(time,Oilcumprod3(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,Oildailyprod3(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% well 4: Oil Production 
figure(7)
title('oil production for well 4')
yyaxis right
plot(time,Oilcumprod4(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,Oildailyprod4(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%


% well 1: water Production 
figure(8)
title('water production for well 1')
yyaxis right
plot(time,watercumprod1(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,waterdailyprod1(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% well 2: water Production 
figure(9)
title('water production for well 2')
yyaxis right
plot(time,watercumprod2(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,waterdailyprod2(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% well 3: water Production 
figure(10)
title('water production for well 3')
yyaxis right
plot(time,watercumprod3(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,waterdailyprod3(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% well 4: water Production 
figure(11)
title('water production for well 4')
yyaxis right
plot(time,watercumprod4(1:n),'r.-');
xlabel('Time(days)','fontsize',13);
ylabel('cumulative production(STB)','fontsize',13);
yyaxis left
plot(time,waterdailyprod4(1:n),'b.-')
xlabel('Time(days)','fontsize',13);
ylabel('Daily production(STB)','fontsize',13);
%

% recovery factor plot
figure(12)
plot(time,R_F(1:n),'r.-')
hand=title('Recovery factor');
set(hand,'fontsize',14);
xlabel('Time(days)','fontsize',13);
ylabel(' Recovery Factor','fontsize',13);


% Water Cut plot
figure(13)
plot(time,f_w(1,1:n),'r.-',time,f_w(2,1:n),'b.-',time,f_w(3,1:n),'g.-',time,f_w(4,1:n),'k.-')
hand=title('Water cut');
set(hand,'fontsize',14);
xlabel('Time(days)','fontsize',13);
ylabel(' Water Cut(%)','fontsize',13);
legend('well 1','well 2','well 3','well 4')

% Average oil phase pressure plot versus time
figure(14)
plot(time,AvPressure(1:n),'r.-')
hand=title('Average reservoir Pressure');
set(hand,'fontsize',14);
xlabel('Time(days)','fontsize',13);
ylabel('Pressure(psi)','fontsize',13);

% Average water saturation  plot versus time
figure(15)
plot(time,AvWaterSat(1:n),'r.-')
hand=title('Average reservoir Water saturation');
set(hand,'fontsize',14);
xlabel('Time(days)','fontsize',13);
ylabel('Water Saturation','fontsize',13);


% presure profiles in some steps
figure(16)
t_break_time=find(f_w(:,:)>0.0001);
t_bt=min(t_break_time);
if t_bt>n
    t_bt=n-1;
end
subplot(2,2,1)
pcolor(POil(:,:,floor(t_bt/4)))
colorbar
hand=title('pressure profile at 1/4 Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

subplot(2,2,2)
pcolor(POil(:,:,floor(2*t_bt/4)));
colorbar
hand=title('pressure profile at 2/4 Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

subplot(2,2,3)
pcolor(POil(:,:,floor(3*t_bt/4)))
colorbar
hand=title('pressure profile at 3/4 Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

subplot(2,2,4)
pcolor(POil(:,:,floor(t_bt)))
colorbar
hand=title('pressure profile at Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

% Water Saturation profiles in some steps
figure(17)
title('dfbrgnbfgn')
subplot(2,2,1)
pcolor(sw(:,:,floor(t_bt/4)))
colorbar
hand=title('Water Saturation profile at 1/4 Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

subplot(2,2,2)
pcolor(sw(:,:,floor(2*t_bt/4)));
colorbar
hand=title('Water Saturation profile at 2/4 Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

subplot(2,2,3)
pcolor(sw(:,:,floor(3*t_bt/4)))
colorbar
hand=title('Water Saturation profile at 3/4 Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv

subplot(2,2,4)
pcolor(sw(:,:,floor(t_bt)))
colorbar
hand=title('Water Saturation profile at Break Through Time');
set(hand,'fontsize',10)
set(gca,'fontsize',10)
colormap hsv