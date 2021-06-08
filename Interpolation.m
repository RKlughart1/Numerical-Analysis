
%% Question 1
clear all;
syms y real
f = @(y)log10(3*y-1);
N = 2;
stepsize = 1/N;
z = 1:stepsize:2;
fx = f(z);

[value,poly]=Lagrange(z,fx,N);%See bottom of file for the function

fprintf("The %d(st/nd/rd/th) order polynomial is: \n",N)
vpa(poly,3)

%% Question 2
clear all;clc;
N = 1:20;

syms y real
f = @(y)log10(3*y-1);

values = [];
RErrors = [];
sigDigs = [];

actual = f(1.4);

for n = 1:20 % find Relative error, significant digits, 
             % and approximation for n = 1:20
    
    stepsize = 1/n;
    z = 1:stepsize:2;
    fx = f(z);
    
    [approx,pol] = Lagrange(z,fx,n);
    values(n) = approx;
    
    relError=abs((actual-approx))/abs(actual);
    RErrors(n) = relError;
    
    sigDig=0;
    while relError<=5*10^(-sigDig)
        sigDig=sigDig+1;
    end
    sigDig=sigDig-1;
    sigDigs(n) = sigDig; 
    
end

table(N',values',RErrors',sigDigs','VariableNames',{'N','Approximation','Relative Error','Significant Digits'})


%% Question 3
%N=5 gives a better approximation than N=19 since one of the values we
%are using to approximate the function is 1.4
clear all
N=5;
x=1:1/N:2;
x(3) % The third element is 1.4

%The interpolating polynomial will be able to predict values better
%if the values we want to predict are close to the values used to
%approximate the polynomial. In this case they are the same value.




%% Lagrange Function

function [value,poly] = Lagrange(z,fx,n)
syms x real;
for i = 1:n+1
    temp=z; % values used for the current iteration
    for j = 1:n+1
        if(j==i)
            target = temp(j);
            temp(j)=[]; % get rid of the element where j==i
        end
    end
    L(i)= prod(x-temp)/prod(target-temp);
end
poly = vpa(sum(fx.*L)); % Interpolating polynomial
value = subs(poly,x,1.4); % evaluating the polynomial at x=1.4

end





