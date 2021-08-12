%% Matlab code for 1D wave equation
% given equation is utt=uxx

clc;clear all;close all;

%% explicit method (intializing the equation)
delx=0.01; %grid resolution size (change this to refine the grid)
r=1; %aspect ratio
u=1; %constant of the wave equation
delt=r^2*delx/2 %time step size
Tsteps= round(1/delt); %noumber of time steps

%% X1 is the potential grid of the simulation due to symmetry only h
%% half of the field is calculated.
X1=zeros(Tsteps,1/(2*delx)+2); %intialize X1

%% Initial condition and reflection line defined
x=0:delx:0.5+delx;
X1(1,:)=sin(pi*x);
X1(2,2:end-1)=0.5*(X1(1,1:end-2)+ X1(1,3:end));
X1(2,end)=X1(2,end-2); %reflection line

for time=3:size(X1,1)
    for space=2:size(X1,2)-1
        X1(time,space)=X1(time-1,space-1)+X1(time-1,space+1)-X1(time-2,space);
    end
    X1(time,end)=X1(time,end-2); %reflected line
end

%% use symmetry condition to create entire field
X2=[X1,fliplr(X1(:,1:end-3))];
AX2=zeros(Tsteps,101); %AX2 is analytical solution
for n=1:1:200
    for x=1:1:101
        AX2(n,x)=sin(pi*x/100)*cos(pi*n/100);
    end
end

%% plotting Numerical solution
figure(1)
imagesc(0:delx:1,(0:delt:Tsteps*delt),X2)
colorbar
ylabel('\leftarrow time(sec)');
xlabel('x');
title('Numerical Hyperbolic PDE')


%% plotting Analytical solution
figure(2)
imagesc(0:delx:1,(0:delt:Tsteps*delt),AX2)
colorbar
ylabel('\leftarrow time(sec)');
xlabel('x');
title('Analytical Hyperbolic PDE')

% %%
% if delx==0.1
%     dispmat=[X1(1:8, 1:7)];
%     fprintf('\n Compare to ta, Solution of wave equation\n');
%     disp(num2str(dispmat))
% end

    
    

