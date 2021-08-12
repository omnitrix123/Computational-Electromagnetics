%% Matlab code to illustrate heat diffusion equation
%% ut=c*uxx first order in time and second order in space)
%% implemented using implicit method
 
clc;close all;clear all

%% Define parameters for diffusion eqn and the range in space and time
L=1.; %length of wire
T=1.; %final time

%% Parameters needed to solve the equation within explicit method
maxk=2500; %no of time steps
dt=T/maxk;
n=50; %no of space steps
nint=50; %The wavefront: intermediate point from which u=0
dx=L/n;
cond=0.25; %conductivity
b=cond*dt/(dx+dx); % stabilty parameter (b<=1)
                      %the stabilty condition depends more on solution
                      %propagating from the point (n,j) hence we use b=2r
                      %for the stability parameter                  
                      
%% Initial temperature of the wire a sinusoid
for i=1:(n+1)
    x(i)=(i-1)*dx;
    u(i,1)=sin(pi*x(i));
end

%% Temperature at the boundary (T=0)
for k=1:maxk+1
    u(1,k)=0.;
    u(n+1,k)=0.;
    time(k)=(k-1)*dt;
end
aa(1:n-2)=-b;
bb(1:n-1)=1.+2.*b;
cc(1:n-2)=-b;
MM=inv(diag(bb,0)+diag(aa,-1)+diag(cc,1));


%% Implementation of the Implicit method
for k=1:maxk %time loop
    for i=2:n %space loop
        u(i,k+1)=u(i,k) +0.5*b*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

%% Graph of temperature at different selected times
figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-');
%at time step 0,100,300,600
title('Temperature within the explicit method');
xlabel('X');
ylabel('T');

figure(2)
mesh(x,time,u');
title('Temperature using the explicit method');
xlabel('X');
ylabel('Temperature');

    