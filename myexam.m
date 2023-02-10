% Trapezoidal rule to find out integral of function 
% f(x) = x^4 - 5*x + 3*(2*x+1)*(x+4)*(7*x+3) - 100 in the interval [3,10]
 
%% initialisation
 
clear all
clc
 
f = @(x) x^4 - 5*x + 3*(2*x+1)*(x+4)*(7*x+3) - 100;  % given function
a = 3;     
b = 10;

g = @(x) (2*x^5 + 105*x^4 + 690*x^3 + 800*x^2 - 640*x)/10 ;   % true integral function of f
trueval = g(10) - g(3)    % calculating true integral value

 
 
n = 100000;
dx = (b-a)/n;
sum = 0;
 
 
%% finding out the value of function for different x value
 
x0 = a;
f0 = f(x0);
 
for i = 1:n
    x(i) = x0 + dx;
    F(i) = f(x(i));
    x0 = x(i);
end
 
%% calculating the integral value
for i = 1:n
    if ( i == 1 || i == n)
        sum = sum + F(i);
    else
        sum = sum + F(i).*2; 
    end
end
It = (dx/2)*sum;
err1 = abs(It - trueval);
disp(['Integral of the function using Trapizoidal rule = ', num2str(It)]);
disp(['Error of integral using Trapizoidal rule = ', num2str(err1)]);


