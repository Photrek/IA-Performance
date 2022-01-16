clc
clearvars
rand('state',0)
randn('state',0)
n=1000;
a=0.0001;b=0.4;B1=2;D1=0.4;D2=1;B2=1;
% I took various values for the degree of freedom, which is equal to the real values of the shape parameter
A=[0.5 0.8 1 1.2 1.5 1.8];

for ia1=1:length(A);
    a1=A(ia1);
    %iteration
for it=1:1000;
    X= trnd(a1,n,1);
    
    Qn=sort(X);
    nn=floor(n*b)-ceil(n*a)+1; 
    
       for jj=1:nn;
    s1(jj)=ceil(n*a)+(jj-1);
end

s1;
s=s1/n;

for j=1:nn
        for ir=1:n
            if (s(j))>((n-ir)/n)&&(s(j))<=((n+1-ir)/n);
            Qt(j)=Qn(ir);
            end
        
        end
    end  
Qt;
Y=log(Qt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0=ones(length(s),1);
X1=-(log(s))';
X2=(2*cos(2*pi*s))';
M=[X1 X0 X2];
W1=s;
W=diag(W1);
M2=inv(M'*W*M)*M'*W*Y';
e=[1 0 0];
Awei=e*M2;
B0_h2(it,1)=Awei;
end
% Finding the mean of the kappa values (because we iteration the experiment (1000) once) 
B0_M2(ia1,1)=mean(B0_h2)
end
Alpha=[0.5 0.8 1 1.2 1.5 1.8]';
No=[1:length(A)]'
Mean2=[B0_M2];
T=table(No,Alpha,Mean2)
