clear all
close all
 clc
format long 
%N=100
 
   M=2 ; N=6
 A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];  % stem

N=length(A);
I=eye(N);

I1=eye(M);



%B=I(:,1:M)
%B=[-0.666549601221895,0.258757856602330;-0.529674783868792,-0.240268472171118;-0.300205473874438,0.714625909757013;-0.430154537230872,-0.603842067126071]



 x0=ones(N,1);
% x0=[1,2,3,4]'
tf=1


% WB0=zeros(N,N);
% ot=0.025;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
% end
%  C0=pinv(WB0); 
 
 

syms t s real
p=simple(expm(A*t));

pf=expm(A*tf);
q=simple(expm(A*(t-s)));
qf=simple(expm(A*(tf-s)));
qft=simple(expm(A*(tf-t)));
xf=pf*x0;



alpha=0.6

iterations=10000
Cost=zeros(iterations,1);
for i=1:iterations

B=rand(N,M);
B=orth(B);

WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 
us=-B'*qf'*C*xf;
u=subs(us,'s',t);
x=simple(p*x0+int(q*B*us, s, 0, t));

Cost(i)=(1-alpha)*int(x'*x,0,tf)+alpha*int(u'*u,0,tf);
i
end 



