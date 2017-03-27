clear all
 clc
 close all
% mp.Digits(32)
% format longG


% load('A.mat')
% A=random_graph(10,0.4);
% N=size(A,1);
M=3;
N=9;
 A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
 A=A+A';
 
% A(1,N)=1;
 

  %A=csvread('E:\network\Electronic networks\s208_st.csvadjacencymatrix.csv');
  
%  A=csvread('E:\network\Food webs\Rhode.csvadjacencymatrix.csv');

%  A = zeros(N,N);



sum1=zeros(N,M); 
sum2=zeros(N,1); 
sum3=zeros(N,1);
tf=1
N=length(A);

Cost11=zeros(M,N)';
iterations=1
for ii=1:iterations

B=rands(M,N)';
MM=M;


B=(sqrt(1*(MM)/trace(B'*B)))*B;

%Xf=expm(A*tf)*expm(A'*tf);

M=min(size((B)));
 

ot=0.01;
Iterations_number=100
Cost=zeros(Iterations_number,1);
Cost3=zeros(Iterations_number,1);
Cost5=zeros(Iterations_number,1);


v=0.03
 
for k=1:Iterations_number
 WB0=zeros(N,N);
 F0=zeros(N,M); 
for k1=1:tf/ot
    WB0=WB0+expm(-A*(ot*k1))*B*B'*expm(-A'*(ot*k1))*ot;
end
 C=pinv(WB0); 
 
 Cost(k)=trace(C);
 
for k1=1:tf/ot
 F0=F0-2*(expm(-A'*(ot*k1))*C*C*expm(-A*(ot*k1)))*B*ot;
end






F1=2*B; % gradient of norm function
F00=F0/norm(F0);  %unit vector of the gradient of energy function 
F11=F1/norm(F1); % unit vector of the gradient of norm function
F11=(sqrt(1*(MM)/trace(F11'*F11)))*F11;    %Normailize gradient of norm function

 
Cost5(k)=sum(F00(:).*F11(:))/(norm(F00(:))*norm(F11(:)));

if  Cost5(k)+1<0.00000005
    break
end

F0=reshape((eye(N*M)-F1(:)*pinv(F1(:)))*F0(:),N,M);

B=B-v*F0/norm(F0);  % F0 is the projected direction 
B=(sqrt(1*(MM)/trace(B'*B)))*B; %initiate the norm to MM


k
   
end

ii
%B=(sqrt(MM))*B/norm(B);

sum1=sum1+(abs(B))/iterations;  % 累加
BB=sum((abs(B)),2);             % 本次

[YY,pos]=sort(BB,'descend');
nodes_leader(:,ii)=pos(1:M);
end

for kk=1:N
e=find(nodes_leader==kk); 
number_leader(kk,1)=length(e);
end 

%% calculate the correlation 
 im=sum(abs(sum1'))/max(sum(abs(sum1')));
 prob=number_leader/(iterations);
 rou= im*prob/(norm(im)*norm(prob))
  
figure(1)

plot(log10(Cost),'r-*')

 
figure(2)
plot(Cost5,'r-*')
% 
% figure(3)
% plot(sum(sum1,2)/max(sum(sum1,2)),'r-*','LineWidth',1.5)% everage of B guiyihua
% 
% 
%  figure(4)
%  plot(number_leader/(iterations),'b-*','LineWidth',1.5) % probability


%% 下面求PGME 按照Science paper的方法求Bs
[N,M]=size(B);
Bb = zeros(N,M);
sum_B=sum(abs(sum1),2);
[max_v order]=sort(sum_B,'descend');%max_v为值，order为在sum_B中的序号
for k=1:M
    Bb(order(k),k)=1;
end

WB0=zeros(N,N);

%先求出文章中的WB，即为这里的WB0
for k1=1:tf/ot
    WB0=WB0+expm(-A*(ot*k1))*Bb*Bb'*expm(-A'*(ot*k1))*ot;
end
C=pinv(WB0);
PGME=trace(C);
log10(PGME)

% rank(ctrb_rank)

rank_control=N-rank(ctrb(A,Bb), 0) % 满zhi为0





%% Random selection method
a=setdiff(1:N,Driver_Nodes);
ND=size(Driver_Nodes,2);
b=randperm(length(a),M-ND);
d0=a(b);
d=[d0,Driver_Nodes];
Bd = zeros(N,M);
for k=1:length(d)
    Bd(d(k),k)=1;
end

WB1=zeros(N,N);
for k1=1:tf/ot
    WB1=WB1+expm(-A*(ot*k1))*Bd*Bd'*expm(-A'*(ot*k1))*ot;
end
C1=inv(WB1);
PGME1=trace(C1);
log10(PGME1)
%pinv
C1=pinv(WB1);
PGME1=trace(C1);
log10(PGME1)
%rank_control=N-rank(ctrb(A,Bd), 0)
u=2*sum(sum(A))/N
