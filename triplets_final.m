clc
clear all
rand('state',0)
randn('state',0)
warning off;
%k = 1.45;
sig=1;
n = 10000;
k=[0.2:0.3:2];
for ia1=1:length(k);
    a1=k(ia1);
for it=1:100
y3=gprnd(a1,1,0,n,1);
A = [];
A1=[];
%Determine triplets from the generated data
for ii=1:1
    r1 = randsample(n, n);
    r2 = randsample(n, n);
    r3 = randsample(n, n);
    A = [A; [y3(r1), y3(r2), y3(r3)]];
end
%Choosing triplets
D=max(A, [], 2)-min(A, [], 2);
[minD, ix] = sort(D);
SS = A(ix,:);
%Choosing (0.5%) from data
S=SS(1:50,:);
% take median
x=median(S,2);
nn=length(x);
%calculating second moment
m=0.0;
b=mean(x);
for i=1:nn
    m=m+(x(i))^2;
end
M_2=(1/nn)*m;
% calculating kappa hat
k_hat=(2*(sig)^2-9*M_2)/3*M_2;%IA papers
hat(it,1)=k_hat;
MS(it,1)=(k_hat-a1)^2;
MML(it,[1 2])=gpfit(x);
end
MSE1(ia1,1)=mean(MS);
B(ia1,1)=mean(hat);
ml(ia1,[1 2])=mean(MML);
va(ia1,1)=var(hat);
end

k=[0.2:0.3:2]';
No=[1:length(k)]';
MSE=[MSE1];Mean1=[B];variance=[va];MLL=[ml];
T=table(No,k,Mean1,variance,MSE,MLL)
