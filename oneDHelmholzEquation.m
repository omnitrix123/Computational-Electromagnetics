%% Matlab code to implement 1D heimholtz equation by FDM
%% It falls under the category of Eigenvalue problem
%% Function 
function k=oneDHelmholzEquation(a,N)

%Arguments:
%a=length of interval
%N=No of subintervals(equal lengths)
%Returns
%k=eigenvalues

h=a/N;%Grid Size
A=spalloc(N-1,N-1,3*(N-1))%allocate sparse matrix with 3*(N-1) zeros
d=-2/h^2; %value of diagonal entries
s=1/h^2; %Value of upper and lower diagonal entries

%% Intialize the diagonal entries
for i=1:N-1
    A(i,i)=d;
end

%% Intialize the upper and lower diagonal entries
for i=1:N-2
    A(i,i+1)=s; %upper diagonal entries
    A(i+1,i)=s; %lower diaonal entries
end

%% computing eigen values 
lambda=eig(A);
k=sqrt(sort(-lambda)); %ort into order of increasing values


    
end


