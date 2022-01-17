clc
clear all
rand('state',0)
randn('state',0)
n=700;
eb=0.001;
a=0.0001;b=0.4;
Beta=2;t=1;
V=[2 1.833 1.667];
for iv=1:length(V);
    v=V(iv);
for it=1:1;
    y=unifrnd(0,1,n,1);
    X = tinv(y,v);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qn=sort(X);
    i=(0:n);
    t=eb+(i/n)*(1-(2*eb));
   for j=1:n+1
        for ir=1:n
        if t(j)>((ir-1)/n)&& t(j)<=(ir/n);
            Qt(j)=Qn(ir);
        end
        end
    end
    
Qt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To build a polynomial
k=n;
for jj=1:(floor(n*b)-ceil(n*a)+1);
    u1(jj)=ceil(n*a)+(jj-1);
end

u1;
u=u1/n;

Le=(1-(2*eb));
Lek=(1/Le^k);

for i1=1:(floor(n*b)-ceil(n*a)+1)
    p1=0.0;
    
    for j1=0:k-1;
 % polynomial function       
        p1=p1+(((Qt(j1+2)-Qt(j1+1))/(1/k))*nchoosek(k-1,j1)*...
                    (u(i1)-eb)^j1*(1-eb-u(i1))^((k-1)-j1));

end
 p(i1)=Lek*p1;
end
p;
Y=-log(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0=ones(length(u),1);
X1=(log(u))';
X2=(2*cos(2*pi*(u)))';
X3=(2*cos(4*pi*u))';
M=[X1 X0 X2 X3];

W1=u/300;
W=diag(W1);
M1=inv(M'*W*M)*M'*W*Y';
e=[1 0 0 0];
Vwei=e*M1;
B0_h2(it,1)=Vwei;
end

B0_M2(iv,1)=mean(B0_h2)
end
V=[2 1.833 1.667]';
No=[1:3]';
Mean2=[B0_M2];
T=table(No,V,Mean2)