%% Ryan Klughart 36627875 A5
clear all; clc;
fprintf("Assignment 5\nRyan Klughart\n36627875\n");
fprintf("\nQuestion 6b\n")
h=0.1;
a=1;
b=2;
fty = @(t,y)y^2/(1+t);
y = [];
t = [];
y(1)=-log(2)^-1;
t(1) = a;

for k = 1:(b-a)/h
    t(k+1) = t(k)+h;
    y(k+1) = y(k) + fty(t(k),y(k))*h;
    fprintf("The approximation at t = %f, is y_%d = %f\n",t(k+1),k,y(k+1));
end
Predicted_b = y;
%%
fprintf("\nQuestion 6d\n")
h=0.1;
a=0;
b=1;
fty = @(t,y)(-t*y)+4*t*y^-1;
y = [];
t = [];
y(1)=1;
t(1) = a;

for k = 1:(b-a)/h
    t(k+1) = t(k)+h;
    y(k+1) = y(k) + fty(t(k),y(k))*h;
    fprintf("The approximation at t = %f, is y_%d = %f\n",t(k+1),k,y(k+1));
end
Predicted_d = y;
%%

fprintf("\nQuestion 8b\n")
actual_y = @(t)(-1)./log(t+1);
a = 1;b=2; h = 0.1; t = a:h:b;
actualValues = actual_y(t);
AE = abs(actualValues - Predicted_b);
table(Predicted_b',actualValues',AE','VariableNames',{'Predicted','Actual','Absolute Error'})

%%
fprintf("\nQuestion 8d\n")

actual_y = @(t)sqrt(4-3*exp(-t.^2));
a = 0;b=1; h = 0.1; t = a:h:b;
actualValues = actual_y(t);
AE = abs(actualValues - Predicted_d);
table(Predicted_d',actualValues',AE','VariableNames',{'Predicted','Actual','Absolute Error'})

%%

fprintf("\nQuestion 14b\n")

fty = @(t,y)(1+t)/(1+y);
a = 1; b = 2; hs = [.5,.1,.02];
y=[];
y(1) = 2;
t=[];
t(1)=a;
for i = 1:3
    fprintf("\nFor h = %f\n",hs(i));
    h= hs(i);
for k = 1:(b-a)/h
    F1 = fty(t(k),y(k));
    ya1 = y(k) + (1/2)*F1*h;
    F2 = fty(t(k)+.5*h,ya1);
    ya2 = y(k) + .5*F2*h;
    F3 = fty(t(k)+.5*h,ya2);
    ya3 = y(k)+F3*h;
    F4 = fty(t(k)+h,ya3);
    y(k+1) = y(k) + (h/6)*(F1+2*F2+2*F3+F4);
    t(k+1) = t(k) + h;
    fprintf("The approximation at t = %f, is y_%d = %f\n",t(k+1),k,y(k+1));
end
end
    
1-800-387-2524.
