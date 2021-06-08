%% Ryan Klughart Assignment 4
%% Question 27
clear all;clc;
syms y real
f = @(y) y^2 + 5*y;
x= 7;
estimations = [];
CDestimations = [];
CD4estimations = [];
for n = 1:20
    h = 10^(-n);
    
    estimate = (f(x+10^(-n))-f(x))/10^(-n);
    estimations(n)=estimate;
    
    CDestimate = (f(x+h)-f(x-h))/(2*h);
    CDestimations(n)=CDestimate;
    
    CD4estimate = (-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);
    CD4estimations(n)=CD4estimate;
    
end

table((1:20)',estimations',CDestimations',CD4estimations','VariableNames',{'N','textbook approximation','CD','CD4'})
% at iteration 12 there is round off error which causes the estimate to
% increase
% at iteration 16, 10^(-16) is rounded to 0 which gives an estimate of 0
% all 3 types of approximations have the same errors at the same iteration

%% Question 7


clear all;clc;
syms y real;
f = @(y) y.^2.*log(y.^2+1);
h = 1/4;
a=0;
b=2;
N = (b-a)/h;%8
x = a:h:b;
fx = f(x);

%Simpsons
Simpsons=0;
for n = 1:N+1
    if(n==1 || n==N+1)
        Simpsons = Simpsons+fx(n);
    end
    if(mod(n,2)==0)
        Simpsons = Simpsons+4*fx(n);
    end
    if(mod(n,2)==1 && n>1 && n<N+1)
        Simpsons = Simpsons+2*fx(n);
    end
end

SimpsonsApprox = Simpsons*(h/3);
Actual = vpa(int(y.^2.*log(y.^2+1),a,b));
    
%Trapezoidal
Trap = 0;
Trap = Trap + f(a)/2 + f(b)/2;
for n = 1:N-1
    Trap = Trap + f(a+n*h);
end
TrapApprox = h*Trap ;

%Midpoint

x = a+h/2:h:b-h/2;
MidpointApprox = h*(sum(f(x)));


Actual = vpa(int(y.^2.*log(y.^2+1),a,b))
TrapApprox = h*Trap
SimpsonsApprox = Simpsons*(h/3)
MidpointApprox = h*(sum(f(x)))


%% Question 14
%% a) Trapezoidal Approximation
clear all;
clc;
syms x real
f = @(x) x.*log(x);
a=1;b=2;
Actual = vpa(int(x*log(x),a,b));
N = 1;
h = (b-a)/N;
x = a:h:b;


Trap = 0;
Trap = Trap + f(a)/2 + f(b)/2;
for n = 1:N-1
    Trap = Trap + f(a+n*h);
end
TrapApprox = h*Trap ;
Ns = [];
hs = [];
errors = [];


Ns(N)=N;
hs(N)=h;
errors(N) = abs(TrapApprox-Actual);
while(abs(TrapApprox-Actual)>10^(-5))
    if(N>1000)
        break;
    end
  
    N=N+1;
    h=(b-a)/N;
    Trap=0;
    Trap = Trap + f(a)/2 + f(b)/2;
    for n = 1:N-1
        Trap = Trap + f(a+n*h);
    end
    TrapApprox = h*Trap ;
    Ns(N)=N;
    hs(N)=h;
    errors(N) = abs(TrapApprox-Actual);
end
fprintf("Table for Trapezoid approximations")
table(Ns',hs',errors','VariableNames',{'N','h','Error'})
fprintf("N = 77 and h = 0.012987 is required for an approximation within 10^-5 for Trapezoid\n")
%% b) Simpsons Approximation
clear all;
clc;
syms x real
f = @(x) x.*log(x);
a=1;b=2;
Actual = vpa(int(x*log(x),a,b));
N = 1;
h = (b-a)/N;
x = a:h:b;
fx = f(x);

Simpsons=0;
for n = 1:N+1 %Find initial approximation
    if(n==1 || n==N+1)
        Simpsons = Simpsons+fx(n);
    end
    if(mod(n,2)==0)
        Simpsons = Simpsons+4*fx(n);
    end
    if(mod(n,2)==1 && n>1 && n<N+1)
        Simpsons = Simpsons+2*fx(n);
    end
end
SimpsonsApprox = Simpsons*(h/3);


Ns = [];
hs = [];
errors = [];
Ns(N)=N;
hs(N)=h;
errors(N) = abs(SimpsonsApprox-Actual);

while(abs(SimpsonsApprox-Actual)>10^(-5))
    if(N>1000)
        break;
    end
    N=N+1;
    h=(b-a)/N;
    x = a:h:b;
    fx=f(x);
    
Simpsons=0;
for n = 1:N+1 
    if(n==1 || n==N+1)
        Simpsons = Simpsons+fx(n);
    end
    if(mod(n,2)==0)
        Simpsons = Simpsons+4*fx(n);
    end
    if(mod(n,2)==1 && n>1 && n<N+1)
        Simpsons = Simpsons+2*fx(n);
    end
end
SimpsonsApprox = Simpsons*(h/3);
Ns(N)=N;
hs(N)=h;
errors(N) = abs(SimpsonsApprox-Actual);
end
fprintf("Table for Simpson approximations")
table(Ns',hs',errors','VariableNames',{'N','h','Error'})
fprintf("N = 6 and h = 0.16667 is required for an approximation within 10^-5 for simpson\n");

%% c) Midpoint Approximation

clear all;
clc;
syms x real
f = @(x) x.*log(x);
a=1;b=2;
Actual = vpa(int(x*log(x),a,b));
N=1;
h = (b-a)/N;
x = a+h/2:h:b-h/2;
MidpointApprox = h*(sum(f(x)));
Ns = [];
hs = [];
errors = [];


Ns(N)=N;
hs(N)=h;
errors(N) = abs(MidpointApprox-Actual);
while(abs(MidpointApprox-Actual)>10^(-5))
    if(N>1000)
        break;
    end
  
    N=N+1;
    h=(b-a)/N;
    x = a+h/2:h:b-h/2;
    MidpointApprox = h*(sum(f(x)));
    
    Ns(N)=N;
    hs(N)=h;
    errors(N) = abs(MidpointApprox-Actual);
end
fprintf("Table for Midpoint approximations")
table(Ns',hs',errors','VariableNames',{'N','h','Error'})
fprintf("N = 54 and h = 0.018519 is required for an approximation within 10^-5 for midpoint\n")
%% 4.6 Question 6c
clear all;clc;
syms y real
Actual = vpa(int(y.*sin(4.*y),-1,1))
a = -1;b = 1;
epsilon = 10^(-5);
AQ(a,b,epsilon)

%% Functions

function approx = AQ(a,b,epsilon)

N=1;
S1= SimpsonsF(a,b,N);
N=2;
S2 = SimpsonsF(a,b,N);

if abs(S2-S1)<epsilon || epsilon<10^-7 % epsilon<10^-7 in case epsilon
    approx=S2;                         % decreases faster than abs(S2-S1)
else
    c = (a+b)/2;
    approx = AQ(a,c,epsilon/2)+AQ(c,b,epsilon/2);
end

end


function approx = SimpsonsF(a,b,N) % to help calculate AQ
syms y real
f = @(y) y.*sin(4.*y);
h=(b-a)/N;
x = a:h:b;
fx=f(x);


approx=0;
for n = 1:N+1 
    if(n==1 || n==N+1)
        approx = approx+fx(n);
    end
    if(mod(n,2)==0)
        approx = approx+4*fx(n);
    end
    if(mod(n,2)==1 && n>1 && n<N+1)
        approx = approx+2*fx(n);
    end
end
approx = approx*(h/3);
end
