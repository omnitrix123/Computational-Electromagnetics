%% Matlab code to solve 2D Poissons equation by FDM
%Numerical scheme used is 5 point differencing method also called as 5
%point stencil

clc;clear all;close all;

%% Setting parameters

lx=10; %length in x direction
ly=10; %length in y direction
nx=10; %no of steps in x
ny=10; %no of steps in y
niter=1000; %no of iterations
dx=lx/(nx-1); %width of space step(x)
dy=ly/(ny-1); %width of space step(y)
x=0:dx:lx; %Range of (0,10) 
y=0:dy:ly;

%% Source functions(as poisson has source present)
f=zeros(nx,ny);
f(round(nx/4),round(ny/4))=3000;
f(round(3*nx/4),round(3*ny/4))=-3000;


%% Intial conditions
p=zeros(ny,nx);
pn=zeros(ny,nx);

%% Boundary conditions(Dirichlet conditions)
p(:,1)=0;
p(:,nx)=0;
p(1,:)=0;
p(ny,:)=0;


%% 5 point stencil
i=2:nx-1;
j=2:ny-1
for it=1:niter
    pn=p;
    p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1)))-(f(i,j)*dx^2*dy^2))/(2*(dx^2+dy^2));
end
    
%% Plotting the solution
surf(x,y,p,'EdgeColor','none');
shading interp
title({'2D Laplace''s equation';['{\itNumber of iterations}=',num2str(it)]})
xlabel('spatial coordinate (x) \rightarrow')
ylabel('{\leftarrow} spatial coordinate (Y)')
zlabel('Solution profile (P) \rightarrow')

    
