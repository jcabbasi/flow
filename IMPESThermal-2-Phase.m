%--- In the Name of God ----
% ---- IMPES MODEL -------
% ------ May,2015 ------------
%--------- BY: Jassem Abbasi and Babak Raji -----------%

clear
clc
%--------- Inputing Data ------------

gn=input('Please Enter Grid Numbers in X direction:');%Grid Numbers
Dx=input('Please Enter Grid Size in X direction:');
Dy=input('Please Enter Grid Size in Y direction:');
Dz=input('Please Enter Grid Size in Z direction:');
dx=Dx*ones(gn);
dy=Dy*ones(gn);
dz=Dz*ones(gn);
Pwp=10000;%BHP of Prod Well 
Pwi=100000;%BHP of Inje Well
muo=2;% Oil Viscosity
muw=1;% Water Viscosity
swinit=0.1;
perm=100; %md
loop=0; % Iteration Counter
ts=1;  
dt=10;
simtsteps=10 %
ministeps=0 %
while ministeps<simtsteps+1
ministeps=ministeps+1;    
re=0.14*((dx(1)^2+dy(1)^2)^0.5);
rw=0.1;
bo=1;
s=1; %skin
k=perm*rand(gn,1);
kro=0.5*ones(gn,1);
so=0.75*ones(simtsteps+1,gn);
loop=0
errfun=1;
%----- Well location selection ---------

wl1=1; %well location
wl2=gn;

%************************************** ITERATION*************************

%-----------------------------------------------------------------------
% errfun: Error Function
while errfun>0.5
% Iterates until reaches Minimum Error
iter=1;
loop=loop+1;
%loop: Number of Iterations
while iter<gn+1
Sw(ts,iter)=1-so(ts,iter);
iter=iter+1;
% Finds Water Saturation for each grid assuming 2 phase flow
end

iter=1;
while iter<gn+1;
sw=Sw(iter);  
%if sw less than critical water saturation:
    if sw<0.15
        krw(iter)=0;
        kro(iter)=0.75;
    elseif sw<0.95
        
krw(iter)= 5.3613*(sw).^4 - 9.3926*(sw).^3 + 4.1554*(sw).^2 + 0.7648*(sw) - 0.1194;
kro(iter)= 9.1783*(sw).^4 - 16.059*(sw).^3 + 7.7899*(sw).^2 - 1.6346*sw + 0.85;  
% if sw more than maximum water saturation
    else
        krw(iter)=0.6;
        kro(iter)=0;
end
iter=iter+1;
end

%---- Well Connection Factor ------------------------
jss1=(((2*pi*k(wl1)*kro(wl1)))/(muo*bo))*(1/(log(re/rw)+s));
jss2=-(((2*pi*k(wl2)*kro(wl2)))/(muo*bo))*(1/(log(re/rw)+s));

%----- Transmisibility Calculation -------------
%landa=landa total
% Here assumes each landa refers to transmissibility beetween each grid and
% ... it's pevious one

iter=1;
landa(1)=0;
while iter<gn
    kavg(iter)=0.5*((1/k(iter+1))+(1/k(iter)));
    
landa(iter+1)= (kavg(iter).*kro(iter))/muo + (kavg(iter).*krw(iter))/muw ;
%sw=swinit*ones(gn,1);
iter=iter+1;
end
landa(gn+1)=0;
%---Clculating Transmissibilities -----
A=zeros(gn+2,gn+2);
iter=1;
while iter<(gn+1)
    A(iter+1,iter+1)=-((landa(iter)/(dx(iter)^2))+(landa(iter+1)/(dx(iter)^2)));
    A(iter+1,iter+1+1)=(landa(iter+1)/(dx(iter)^2));
    A(iter+1,iter-1+1)=landa(iter)/(dx(iter)^2);
    
    
    iter=iter+1;
end


iter=1;
while iter<gn+1
    AA(iter,iter)=A(iter+1,iter+1);
   iter=iter+1; 
end
iter=1;
while iter<gn
    AA(iter+1,iter)=A(iter+2,iter+1);
   iter=iter+1; 
end

iter=1;
while iter<gn
    AA(iter,iter+1)=A(iter+1,iter+2);
   iter=iter+1; 
end

B=zeros(gn,1);
%---- Well Effect in Pressure Decline:
AA(wl1,wl1)=AA(wl1,wl1)- jss1;
B(wl1,1)= jss1*Pwp;

AA(wl2,wl2)=AA(wl2,wl2)- jss2;
B(wl2,1)= jss2*Pwi;
% ----------------- Finding Pressure at each grid:
X=inv(AA)*B;
X=abs(X);
%----------------------Saturation Calculation ---------------

while iter<gn+1
Sw(ts,iter)=1-so(ts,iter);
iter=iter+1;

end




iter=1;
while iter<gn   
%kavg(iter)=0.5*((1/k(iter+1))+(1/k(iter)));%Recalculation    
landao(iter+1)= (kavg(iter).*kro(iter))/muo;
%sw=swinit*ones(gn,1);% I guess it is a problem source
iter=iter+1;
end

alfa=1;
iter=1;

while iter<gn+1
    if iter==1
 so(ts+1,iter)= so(ts,iter)+ (dt/alfa)*(((landao(iter)*(X(iter+1,1)-X(iter,1)))/dx(iter)^2));
    elseif  iter==gn
 so(ts+1,iter)= so(ts,iter)- (dt/alfa)*(((landao(iter-1)*(X(iter)-X(iter-1)))/dx(iter)^2) );
    else
 so(ts+1,iter)= so(ts,iter)+ (dt/alfa)*(((landao(iter)*(X(iter+1,1)-X(iter,1)))/dx(iter)^2)- ((landao(iter-1)*(X(iter)-X(iter-1)))/dx(iter)^2) );
    end     
iter=iter+1;
end
iter=1;
while iter<gn+1
errsat(iter)= abs(so(ts+1,iter)-so(ts,iter));
iter=iter+1;
end
iter=1;
while iter<gn+1
  if errsat(iter)>0.000000000005
      errfun=errsat(iter);
  iter=iter+1;
  else
      errfun=0;
iter=iter+1;
  end
  
end
end
Iteratenumbers(ministeps,1)=loop
ts=ts+1
end

loop;