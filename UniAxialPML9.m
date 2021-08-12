%% Matlab code to study scattering of EM waves from an
%% object where the domain is truncated using UPML
clc;close all;clear all;

%% Physical domain
a=50;
b=50;

%% Constants
eps0=8.852e-12; %permittivity of free space
mu0=4*pi*1e-7; %permeability of free space

%% Source parameters
c=3e8; %Speed of EM wave
freq=3e7; %frequency of EM wave=300 MHz
lambda=c/freq;
k0=(2*pi)/lambda;
w=2*pi*freq;

%% Grid Parameters
dx=lambda/20; %mesh size along X direction
dy=lambda/20; %mesh size along Y direction

x=0:dx:a;
nx=length(x);

y=0:dy:b;
ny=length(y);

%% Source location assigned in middle
%floor(nx/2),floor(ny/2)

%% Time parameters
sim_time=200;
R=0.5;
dt=R*dx/c;
tsteps=sim_time;

%% Compute PML parameters
pmlp=2; %profile of conductivity inside PML
refl_th=-70 % Theoretical numerical reflection dB

%sig0x=((eps0*c)/(d_pmlx))*((pmlprof+1)/2)*ln(1/reftheory);
sig0x=0.01;
sig0y=0.01;

d_xpml=8;
d_ypml=8;
leftrightpml=[1 1];
lowuppml=[1 1];

sigx=zeros(nx,ny);
sigy=zeros(nx,ny);
mu=zeros(nx,ny);
eps=zeros(nx,ny);
z=zeros(nx,ny);

mur=ones(nx,ny);
epsr=ones(nx,ny);
oneM=ones(nx,ny);
dtM=dt*oneM;

%% Filling up the conductivity values
for i=1:nx
    for j=1:ny
        if (leftrightpml(1)>0) && (i<= d_xpml)
            sigx(i,j)=sig0x*(((d_xpml-i)*dx)/(d_xpml*dx))^pmlp;
        end
        
        if (leftrightpml(2)>0) && (i>=(nx-d_xpml))
            sigx(i,j)=sig0x*(((i-(nx-d_xpml))*dx)/(d_xpml*dx))^pmlp;
        end
        
        if (lowuppml(1)>0) && (j<=d_ypml)
            sig(i,j)=sig0y*(((d_ypml-j)*dy)/(d_ypml*dy))^pmlp;
        end
        
        if (lowuppml(2)>0) && (j>=(ny-d_ypml))
            sig(i,j)=sig0y*(((j-(ny-d_ypml))*dy)/(d_ypml*dy))^pmlp;
        end 
        
        mu(i,j)=mur(i,j)*mu0;
        eps(i,j)=epsr(i,j)*eps0;
        z(i,j)=sqrt(mu(i,j)/eps(i,j));
    end
end


%% Initialize the field matrices
Hx=zeros(nx,ny);
Hy=zeros(nx,ny);
Ez=zeros(nx,ny);
kxx=zeros(nx,ny);
kyy=zeros(nx,ny);

%% Intialize curl matrices
CHx=zeros(nx,ny); 
CHy=zeros(nx,ny);
CEz=zeros(nx,ny);

%% Initialize integration matrices
ICHx=zeros(nx,ny);
ICHy=zeros(nx,ny);
ICEz=zeros(nx,ny);
ICkx=zeros(nx,ny);
ICky=zeros(nx,ny);

%% Starting the main FDTD loop for E field
for t=1:tsteps
    
    %Boundary condition PEC
    Hx(:,1)=0;
    Hx(:,ny)=0;
    Hy(1,:)=0;
    Hy(nx,:)=0;
    
    
    Ez(:,1)=0;
    Ez(:,ny)=0;
    Ez(1,:)=0;
    Ez(nx,:)=0;
    
    
    %Defining a scatterer
    for i=73:74
        for j=48:52
            Ez(i,j)=0;
        end
    end
    
    % Calculating CHx
    for i=1:nx
        for j=1:ny-1
            kxx(i,j)=kxx(i,j)-(dt*((Ez(i,j+1)-Ez(i,j))/dy));
            
            Hx(i,j)=Hx(i,j)-(dt/mu(i,j))*((Ez(i,j+1)-Ez(i,j))/dy)...
                +(((dt*sig(i,j))/(eps(i,j)*mu(i,j)))*kxx(i,j))...
                -(dt*sigy(i,j)/mu(i,j))*Hx(i,j);
            
        end
    end
    
    %Calculating CHy
    for i=1:nx-1
        for j=1:ny
            kyy(i,j)=kyy(i,j) + (dt*((Ez(i+1,j)-Ez(i,j))/dx));
            
            Hy(i,j)=Hy(i,j)+(dt/mu(i,j))*((Ez(i,j))/dx)...
              +(dt/(eps(i,j)*mu(i,j))*sigy(i,j)*kyy(i,j))...
              - (dt/mu(i,j))*sigx(i,j)*Hy(i,j);
         
        end
    end
    
    for i=2:nx
        for j=2:ny
            Ez(i,j)=Ez(i,j)+(dt/eps(i,j))*((Hy(i,j)-Hy(i-1,j))/dx-(Hx(i,j)-Hx(i,j-1))/dy);
        end
    end
    
    Ez(1,1)=Ez(1,1)+(dt/eps(1,1))*((Hy(1,1)+Hy(1,1))/dx-(Hx(1,1)+Hx(1,1))/dy);
    
    %source 
    source=sin(((2*pi*freq*t*dt)));
    
    Ez(floor(nx/2),floor(ny/2))=2*source;
    
    %plotting Ez wave
    [xx,yy]=meshgrid(x,y);
    pcolor(x,y,Ez');
    rectangle('position',[36 23.5 2 2], 'facecolor','m');
    shading interp;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    titlestring=['\fontsize{18} FDTD setup with UMPL at timestep = ',num2str(t)];
    axis([0 a 0 b -1 1]);
    view(0,90);
    caxis([-1,1])
    colormap(jet);
    colorbar;
    getframe;
end

    