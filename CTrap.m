% author: Venkata Ayyalasomayajula

% Composite-Trapezoidal method

         
% f     : integrand
% a     : lower limit of integration
% b     : upper limit of integration
% n     : number of sub-intervals

% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=CTrap(f,a,b,n)

%verify if n is an integer and a,b are real numbers
if isreal(a)==0 
    error('Please enter a real number for the lower limit')
elseif isreal(b)==0
    error('Please enter a real number for the upper limit')
elseif rem(n,1)~=0
    error('Please enter an integer value for the number of sub-intervals')
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

for i=1:n-1
    w(i+1)=2;
end

%function evaluation at nodes
f1=feval(f,x);

%evaluation of integral
I=(h/2)*sum(f1.*w);

%plotting
plot(x,f1,'LineWidth',2)
xlabel('x','FontSize',14)
ylabel('f(x)','FontSize',14)
end
