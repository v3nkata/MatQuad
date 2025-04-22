% author: Venkata Ayyalasomayajula

% Gauss-Chebyshev quadrature (first kind)

% w(x)  : 1/sqrt(1-x^2)
% domain: (-1,1)

% f     : integrand
% n     : number of integration points


% x     : nodes
% w     : weights
% I     : approximate value of the integral

function [I,x,w]=GCheby(f,n)

%verify if n is an integer
if rem(n,1)~=0
    error('Please enter an integer value for the number of integration points')
end

%nodes
x=cos(((1:2:2*n)*pi)/(2*n));     %values are in the order of largest-->smallest
x=x(end:-1:1)';               

%weights
w=zeros(n,1);
for i=1:n
    w(i)=pi/n;
end

%evaluate function at the nodes
f1=feval(f,x);

%evaluate quadrature
I=sum(w'*f1);
