% author: Venkata Ayyalasomayajula

% Gauss-Laguerre quadrature 

% w(x)  : exp(-x)
% domain: [0,inf)

% f     : integrand
% n     : number of integration points


% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=GLagu(f,n)

%verify if n is an integer
if rem(n,1)~=0
    error('Please enter an integer value for the number of integration points')
end

%nodes
a=2*(0:n-1)+1;                          %main diagonal
b=(1:n-1);                              %first diagonal above the main diagonal

J=diag(b,1)+diag(a,0)+diag(b,-1);       %tridiagonal matrix J

[V,l]=eig(J);
[x,i]=sort(diag(l));

%weights
w=V(1,i).^2;

%evaluate function values
f1=feval(f,x);

%evaluate quadrature
I=sum(w'.*f1);
