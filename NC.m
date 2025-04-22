% author: Venkata Ayyalasomayajula

% Newton-Cotes quadrature

% f     : integrand
% a     : lower limit of integration
% b     : upper limit of integration
% n     : number of sub-intervals
% m     : degree of quadrature

% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=NC(f,a,b,n,m)

%verify if a,b,n are integers
if isreal(a)==0 
    error('Please enter a real number for the lower limit')
elseif isreal(b)==0
    error('Please enter a real number for the upper limit')
elseif rem(n,1)~=0
    error('Please enter an integer value for the number of sub-intervals')
elseif rem(m,1)~=0
    error('Please enter an integer value for the degree of quadrature')
end

h=(b-a)/(m*n); %width

%nodes
x=a:h:b;

%procedure to compute weights

%coefficient matrix of the system of equations  
V=fliplr(vander(1:m+1))'; 

%vector of constant terms
M=zeros(1,m+1); %pre-allocation
for i=1:m+1
    M(i)=((m+1)^i-1)/i;
end

%weights
w=(V\M')';

%procedure to compute quadrature
f1=feval(f,x);

f2=zeros(m+1,n); %pre-allocation

for i=1:n
    f2(1:m+1,i)=f1(m*(i-1)+1:i*m+1);
end

I1=h*w*f2;
I=sum(I1);




