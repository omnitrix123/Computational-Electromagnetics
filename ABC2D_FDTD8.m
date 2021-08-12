%% Matlab code to demonstrate 2D absorbing boundary condition
% close all;clear all;clc;

%% Define parameters for the simulation
% computational domain in x (xdim) direction

c=3e8; %spped of light
freq_in=3e7; %frequency in Hz
eps_r=1; %relative dielectric constant of medium

lambda=(c/freq_in)/sqrt(eps_r);
xdim=100;
dx=lambda/10;
x=0:dx:xdim;
xsteps=length(x);


%% Total simulation time
time_tot=400;
Rx=0.5; % courant constant of stability<1;
dt=Rx*dx/c;
tsteps=time_tot;

%% Grid Dimension in y (ydim) direction
ydim=100;
Ry=0.5;
dy=c*dt/Ry;
y=0:dy:ydim;
ysteps=length(y);


%% position of source
xsource=floor(xsteps/2);
ysource=floor(ysteps/2);

%% Initialization of field vectors
Ez=zeros(ysteps,xsteps);
Hx=zeros(ysteps,xsteps);
Hy=zeros(ysteps,xsteps);

Ex2=zeros(tsteps,xsteps);
Exlast_1=zeros(tsteps,xsteps);
Ey2=zeros(tsteps,ysteps);
Eylast_1=zeros(tsteps,ysteps);

for n=1+ceil((1/min(Rx,Ry))):tsteps
    
    for l=1:xsteps 
        for m=1:ysteps-1
            Hx(m,l)=Hx(m,l)-Ry*(Ez(m+1,l)-Ez(m,l));
        end
    end
    
    for m1=1:ysteps 
        for l1=1:xsteps-1
            Hy(m1,l1)=Hy(m1,l1)+Rx*(Ez(m1,l1+1)-Ez(m1,l1));
        end
    end
    
    for m2=2:ysteps 
        for l2=2:xsteps
            Ez(m2,l2)=Ez(m2,l2)+Rx*(Hy(m2,l2)-Hy(m2,l2-1))-Ry*(Hx(m2,l2)-Hx(m2,l2-1));
        end
    end
    
    
    %% Absorbing Boundary Condition- using one way  wave equation
    %% Performance deteriorates for off normal incidence
    
    % in x direction
    
    Ex2(n,:)=Ez(:,2);
    Ez(:,1)=Ex2(n-1/Rx,:);
    Exlast_1(n,:)=Ez(:,xsteps-1);
    Ez(:,xsteps)=Exlast_1(n-1/Rx,:);
    
    % in y direction
    Ey2(n,:)=Ez(2,:);
    Ez(1,:)=Ey2(n-1/Ry,:);
    Eylast_1(n,:)=Ez(ysteps-1,:);
    Ez(ysteps,:)=Eylast_1(n-1/Ry,:);
    
    %% Define Source
    source=sin(((2*pi*(freq_in)*n*dt)));
    
    %% Assigning Source 
    Ez(ysource,xsource)=source;
    
    %% plotting Ez-Wave
    mesh(x,y,Ez,'LineWidth',2);
    xlabel('X');
    ylabel('Y');
    zlabel('Ez');
    titlestring=['\fontsize{20}Plot of Ez vs X and Y for 2D FDTD open boundary at time step=',num2str(n)]
    title(titlestring,'color','k')
    axis([0 100 0 100 -1 1]);
    %colormap(jet);
    colorbar;
    getframe;
end

    