% author: Venkata Ayyalasomayajula

% Composite-Simpson method

% f     : integrand
% a     : lower limit of integration
% b     : upper limit of integration
% n     : number of sub-intervals (MUST be even)

% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=CSimp(f,a,b,n)

%verify if n is an integer and a,b are real numbers
if isreal(a)==0 
    error('Please enter a real number for the lower limit')
elseif isreal(b)==0
    error('Please enter a real number for the upper limit')
elseif rem(n,1)~=0
    error('Please enter an integer value for the number of sub-intervals')
elseif rem(n,2)~=0
    error('Please enter an even number for the number of sub-intervals') 
end

h=(b-a)/(n); %width

%nodes
x=zeros(1,n+1); %pre-allocation, to avoid iterative resizing
x(1)=a; 

for k=1:n
    x(k+1)=a+k*h;
end

%weights
w=zeros(1,n+1); %pre-allocation, to avoid iterative resizing

w(1)=1;
w(n+1)=1;

w(2)=4;
for i=2:n-2
    w(i+2)=4;
end

w(3)=2;
for j=3:n-3
    w(j+2)=2;
end

%function evaluation at nodes
f1=feval(f,x);

%evaluation of integral
I=(h/3)*sum(f1.*w);


end
