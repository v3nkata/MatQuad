% author: Venkata Ayyalasomayajula

% Gauss-Lobatto quadrature 

% w(x)  : 1
% domain: [-1,1]

% f     : integrand
% n     : number of integration points
% a     : lower limit of integration
% b     : upper limit of integration

% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=GLoba(f,a,b,n)

%verify if n is an integer and a,b are real numbers
if rem(n,1)~=0
    error('Please enter an integer value for the number of integration points')
elseif isreal(a)==0
    error('Please enter a real number for the lower limit of integration')
elseif isreal(b)==0
    error('Please enter a real number for the upper limit of integration')
end

%nodes
N=n-2;
b1=sqrt(1./(4-1./(1:N).^2));            % first diagonal above the main diagonal
J1=diag(b1,1)+diag(b1,-1);
J=[J1,[zeros(N,1);sqrt((N+1)/(2*N+1))];[zeros(1,N),sqrt((N+1)/(2*N+1))],0];         %additional terms in the modified Jacobi matrix

[V,l]=eig(J);
[x1,i]=sort(diag(l));

x=(b-a)/2*x1+(a+b)/2;                     % change of interval [a,b]-->[-1,1]

%weights
w1=2*V(1,i).^2;

w=(b-a)/2*w1;                             % change of interval [a,b]-->[-1,1]

%evaluate function value at nodes
f1=feval(f,x);

%evaluate quadrature
I=sum(w'.*f1);
