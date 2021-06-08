clear all;clc;
fprintf("\nQuestion 6c\n")
A = [4 1 -1 1; 1 4 -1 -1; -1 -1 5 1; 1 -1 1 3];
b = [-2 -1 0 1]';
size(b)
E = [A b];
D = diag(A)';
D = diag(D);
M = A-D;
x0 = [0 0 0 0]';
x = [];
x(1,:) = x0;
n = 2;
for k =1:n
    x(k+1,:) = inv(D)*(b-M*x(k,:)');
    fprintf("%d(st/nd) iteration estimate x=",k)
    x(k+1,:)'
end

%% Question 8c
fprintf("\nQuestion 8c\n")
N = 2; n = 4;
x0 = [0 0 0 0]';
x = [];
x(1,:) = x0;
L =tril(A);
R = A-L;
for k = 1:N
    x(k+1,:) = inv(L)*(b-R*x(k,:)');
     fprintf("%d(st/nd) iteration estimate x=",k)
    x(k+1,:)'
end
        
%%

fprintf("\nQuestion 14b\n")

ns = [10,50,100,9,99];
for k = 1:length(ns)
    n=ns(k);
    A = zeros(n,n);
    b = zeros(n);
    b(1) = 1/2;
 for i = 1:n
     for j = 1:n
         if(i==j)
             A(i,j)=1;
             if(i<n)
             A(i+1,j)=-1/2;
             A(i,j+1)=-1/2;
             end
         end
     end
 end

 P = A\b(:,1);
 AP = A*P; % = B
 fprintf("For n = %d\n",n);
table(round(P,4),round(AP,4),'VariableNames',{'P','AP'})
end




%% Question 14d

fprintf("\nQuestion 14d\n")
alpha = 0.3;
ns = [10,50,100,9,99];
for k = 1:length(ns)
    n=ns(k);
    A = zeros(n,n);
    b = zeros(n);
    b(1) = 1/2;
 for i = 1:n
     for j = 1:n
         if(i==j)
             A(i,j)=1;
             if(i<n)
             A(i+1,j)=alpha;
             A(i,j+1)=1-alpha;
             end
         end
     end
 end
 P = A\b(:,1);
 fprintf("For n = %d\n",n);

 AP = A*P; % = B

 table(round(P,4),round(AP,4),'VariableNames',{'P','AP'})
end

%% Question 14 b p2


fprintf("\nQuestion 14b partII\n");


ns = [10,10^2,10^3,10^4];
for k = 1:length(ns)
    tic
    n=ns(k);
    x = [];
    x0 = zeros(n,1);
    x(1,:) = x0;
    A = zeros(n,n);
    b = zeros(n,1);
    b(1) = 1/2;
 for i = 1:n
     for j = 1:n
         if(i==j)
             A(i,j)=2;
             if(i<n)
             A(i+1,j)=-.5;
             A(i,j+1)=-.5;
             end
         end
     end
 end
 
D = diag(A)';
D = diag(D);
M = A-D;


for k =1:n
    x(k+1,:) = D\(b-M*x(k,:)');
    if(norm(x(k+1,:)-x(k,:))<10^-3)
       
        fprintf("Solution found at iteration %d using Jacobi\n",k);
         toc
        break;
    end
end

end



ns = [10,10^2,10^3,10^4];
for k = 1:length(ns)
    tic
    n=ns(k);
    x = [];
    x0 = zeros(n,1);
    x(1,:) = x0;
    A = zeros(n,n);
    b = zeros(n,1);
    b(1) = 1/2;
 for i = 1:n
     for j = 1:n
         if(i==j)
             A(i,j)=2;
             if(i<n)
             A(i+1,j)=-.5;
             A(i,j+1)=-.5;
             end
         end
     end
 end
 
L =tril(A);
R = A-L;
for k = 1:n
    x(k+1,:) = L\(b-R*x(k,:)');
    if(norm(x(k+1,:)-x(k,:))<10^-3)
        fprintf("Solution found at iteration %d using GS\n",k);
        toc
        break;
    end
end


end

%% Question 14
fprintf("\nQuestion 14d partII\n");

alpha=0.3;

%Jacobi
ns = [10,10^2,10^3,10^4];
for k = 1:length(ns)
    tic
    n=ns(k);
    x = [];
    x0 = zeros(n,1);
    x(1,:) = x0;
    A = zeros(n,n);
    b = zeros(n,1);
    b(1) = 1/2;
 for i = 1:n
     for j = 1:n
         if(i==j)
             A(i,j)=2;
             if(i<n)
             A(i+1,j)=alpha;
             A(i,j+1)=1-alpha;
             end
         end
     end
 end
 
D = diag(A)';
D = diag(D);
M = A-D;


for k =1:n
    x(k+1,:) = D\(b-M*x(k,:)');
    if(norm(x(k+1,:)-x(k,:))<10^-3)
       
        fprintf("Solution found at iteration %d using Jacobi\n",k);
         toc
        break;
    end
end

end

%GS

ns = [10,10^2,10^3,10^4];
for k = 1:length(ns)
    tic
    n=ns(k);
    x = [];
    x0 = zeros(n,1);
    x(1,:) = x0;
    A = zeros(n,n);
    b = zeros(n,1);
    b(1) = 1/2;
 for i = 1:n
     for j = 1:n
         if(i==j)
             A(i,j)=2;
             if(i<n)
             A(i+1,j)=alpha;
             A(i,j+1)=1-alpha;
             end
         end
     end
 end
 
L =tril(A);
R = A-L;
for k = 1:n
    x(k+1,:) = L\(b-R*x(k,:)');
    if(norm(x(k+1,:)-x(k,:))<10^-3)
        fprintf("Solution found at iteration %d using GS\n",k);
        toc
        break;
    end
end


end

       
%% Question 6c
fprintf("\nQuestion 6c\n")
A = [4 1 -1 1; 1 4 -1 -1; -1 -1 5 1; 1 -1 1 3];
b = [-2 -1 0 1]';
w = 1.2;
tol = 10^-3;
k= 1;
N=100;
n = length(b);
x = zeros(n,1);
XO = zeros(n,1);
while(k<N)
    for(i = 1:n)
        
        sum1 = 0;
        for(j = 1:i-1)
            if(i~=1)
                sum1 =sum1+ A(i,j).*x(j,:);
            end
        end
        
        sum2 = 0;
        for(j = i+1:n)
            sum2 = sum2+ A(i,j)*XO(j,:);
        end
        
        x(i) = (1-w)*(XO(i)) + 1/(A(i,i))*w*(-sum1-sum2+b(i));
        
        if(norm(x-XO)<tol)
            fprintf("Solution found at iteration %d\n",k);
            x
            k=N;
            break;
        end
    end
    k=k+1;
    for(i = 1:n)
        XO(i) = x(i);
    end
end

