%% Matlab code for 1D electromagnetic FDTD program
%% Assumes Ex and Hy field  components propagating in the
%% x direction. 
%% FIELDS, PERMITTIVITY, PERMEABILITY AND CONDUCTIVITY ARE FUNCTIONS OF X

clc;close all;clear all;

%% Specifying the parameters

L=5; %domain length
N=505; % # of spatial samples in domain
Niter=800; % # of iterations to perform
fs=300e6; %source frequency in Hz

ds=L/N; %spatial step in Hz
dt=ds/300e6; %magic time step


eps0=8.854e-12; %permittivity of free space
mu0=pi*4e-7; %permeability of free space

x=linspace(0,L,N); %x coordinate of spatial samples
showWKB=0; % if =1 then show WKB approx at end

%% Scale factors for E and H
ae=ones(N,1)*dt/(ds*eps0);
am=ones(N,1)*dt/(ds*mu0);
as=ones(N,1);
epsr=ones(N,1);
mur=ones(N,1);
sigma=zeros(N,1);



%% Profile sequence for the epsilion
%profile sequence 1,2 3,4,

profile=2;

%% different profile definitions

for i=1:N
    epsr(i)=1;
    mur(i)=1;
    w1=0.5;
    w2=1.5;
    
    if (profile==1)
        if (abs(x(i)-L/2)<0.5)
            epsr(i)=4;
        end
    end
    
    if (profile==2)
        if(abs(x(i)-L/2)<1.5)
            epsr(i)=1+3*(1+cos(pi*(abs(x(i)-L/2))))
        end
        
        if (abs(x(i)-L/2)<0.5)
             epsr(i)=4;
        end
    end
     
    if(profile==3)
        if (x(i)>L/2) 
            epsr(i)=9;
        end
    end
    
    if(profile==4)
        if(x(i)>(L/2-0.1443))
            epsr(i)=3;
        end
        
        if(x(i)>L/2)
            epsr(i)=9;
        end
    end
    
    if(profile==5)
        if(x(i)>L/2)
            sigma(i)=0.005;
        end
    end
    
    if(profile==6)
        epsr(i)=1+sin(2*pi*x(i)/L)^2;
    end
    
end

ae=ae./epsr;
am=am./mur;
ae=ae./(1+dt*(sigma./epsr)/(2*eps0));
as=(1 -dt*(sigma./epsr)/(2*eps0))./(1+dt*(sigma./epsr)/(2*eps0));

%% Plot the permittivity permeability and conductivity profiles

figure(1)
subplot(311)
plot(x,epsr);
grid on;
axis([3*ds L min(epsr)*0.9 max(epsr)*1.1]);
title('Relative permittivity');

subplot(312)
plot(x,mur);
grid on;
axis([3*ds L min(mur)*0.9 max(mur)*1.1]);
title('Relative permeability');

subplot(313)
plot(x,sigma);
grid on;
axis([3*ds L min(sigma)*0.9-0.001 max(sigma)*1.1+0.001]);
title('Conductivity');


