%% Matlab code for calculating the capacitor of a pair
%% of coaxial rectangles using FDM
%% Governing equations are the 2D Laplaces Equation

function cap=capacitance(a,b,c,d,n,tot)

%Arguments:
% a width of inner conductor
% b height of inner conductor
% c width of outer conductor
% d height of outer conductor
% n no of pointd in the x direction(horizontal)
% tot relative tolerance for capacitance
% rel relaxation parameter=2-c/n

%Returns:
% cap = capacitance per unit length (pf/m)

%%  Make grid
rel=2-c/n;               %Relaxation parameter
h=0.5*c/n;               %Grid size
na=round(0.5*c/h);       %Number of segments on a
x=linspace(0,0.5*c,n+1); %Grid points along x axis
m=round(0.5*d/h);        %Number of segments on d
mb=round(0.5*b/h);       %Number of segments od b
y=linspace(0,0.5*d,m+1); %Grid points along y axis

%% Initialize potential and mask array
 f=zeros(n+1,m+1);       %2D array with solution
 mask=ones(n+1,m+1)*rel; %2D array with relaxation
 
 for i=1:na+1
     for j=1:mb+1
         mask(i,j)=0;
         f(i,j)=1;
     end
 end
 
 
 %% Gauss Seidel Iteration
 oldcap=0;
 for iter =1:1000    %Maximum no of iterations
     f=seidel(f,mask,n,m); %Perform Gauss Seidel Iteration
     cap=gauss(n,m,h,f);  %Compute the Capacitance
     if(abs(cap-oldcap)/cap<tol)
         break %stop if change in capacitance is sufficiently small
     else
         oldcap=cap; %continue untill converged
     end
 end
 
end

        


