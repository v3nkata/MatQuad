% author: Venkata Ayyalasomayajula

% Gauss-Hermite quadrature 

% w(x)  : exp(-x^2)
% domain: (-inf,inf)

% f     : integrand
% n     : number of integration points


% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=GHerm(f,n)

%verify if n is an integer
if rem(n,1)~=0
    error('Please enter an integer value for the number of integration points')
end

%nodes
b=sqrt((1:n-1)/2);              %first diagonal above the main diagonal
J=diag(b,1)+diag(b,-1);         %tridiagonal Jacobi matrix
[V,l]=eig(J);
[x,i]=sort(diag(l));            

%weights
w=sqrt(pi)*V(1,i)'.^2;

%evaluate function at nodes
f1=feval(f,x);

%evaluate quadrature
I=sum(w'*f1);
