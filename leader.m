clear all
 clc
 close all
 

load ('A.mat')
N=size(A)
M=25
%NN=5
%  A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
% A(1,N)=1;

% A good example to show the minimization of  the maximum matching paths
% corresponds to the optimal energy control.
 
%A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
  %A=csvread('E:\network\Electronic networks\s208_st.csvadjacencymatrix.csv');
  
%  A=csvread('E:\network\Food webs\Rhode.csvadjacencymatrix.csv');

% A = zeros(N,N);

% load('D:\code\dataset\Foodweb\Rhode.mat')
% L = size(A,1)%L为边的个数
% N2 = length(unique(A))
% [A,~] = Net_Consecutive(double(A)); 
% N=size(A);
% A=A';



% 
%   for i=1:N-1
%      A(i+1,i)=1;
%   end % stem
 
%A(1,N)=1 %cycle
% 
%  A(NN,1)=1; %dialtion
%  A(NN,NN-1)=0;
%  


sum1=zeros(N,M); 
sum2=zeros(N,1); 
sum3=zeros(N,1);
tf=2
N=length(A);

Cost11=zeros(M,N)';
iterations=1
for ii=1:iterations

B=rands(M,N)';
MM=M;


B=(sqrt(1*(MM)/trace(B'*B)))*B;

Xf=expm(A*tf)*expm(A'*tf);

M=min(size((B)));
  

ot=0.025;

 
Iterations_number=500;
Cost=zeros(Iterations_number,1);
Cost3=zeros(Iterations_number,1);
Cost5=zeros(Iterations_number,1);
v=0.03
 
for k=1:Iterations_number
 WB0=zeros(N,N);
 F0=zeros(N,M); 
for k1=1:tf/ot
    WB0=WB0+expm(A*(ot*k1))*B*B'*expm(A'*(ot*k1))*ot;
end
 C=pinv(WB0); 
for k1=1:tf/ot
 F0=F0-(expm(A'*(ot*k1))*C*Xf*C*expm(A*(ot*k1)))*B*ot;
end

 
Cost(k)=trace(C*Xf);



F1=2*B; % gradient of norm function
F00=F0/norm(F0);  %unit vector of the gradient of energy function 
F11=F1/norm(F1); % unit vector of the gradient of norm function
F11=(sqrt(1*(MM)/trace(F11'*F11)))*F11;    %Normailize gradient of norm function

Cost3(k)=trace(B'*B);
 
Cost5(k)=sum(F00(:).*F11(:))/(norm(F00(:))*norm(F11(:)));

if  Cost5(k)+1<0.000000005
    break
end

F0=reshape((eye(N*M)-F1(:)*pinv(F1(:)))*F0(:),N,M);

B=B-v*F0/norm(F0);  % F0 is the projected direction 
B=(sqrt(1*(MM)/trace(B'*B)))*B; %initiate the norm to MM


k
   
end
BBB_tensor(:,:,ii)=B;
ii
B=(sqrt(MM))*B/norm(B);
sum1=sum1+(abs(B))/iterations;
BB=sum((abs(B)),2);
%sum2=sum2+(BB/iterations)/sum(BB);
sum2=sum2+(BB/iterations)/max(BB);

Cost11=Cost11+B/iterations;

[YY,pos]=sort(BB,'descend'); % this iteration selected leader
nodes_leader(:,ii)=pos(1:M);
end

for kk=1:N
e=find(nodes_leader==kk); 
number_leader(kk,1)=length(e);
end 

figure(1)

plot(Cost,'r-*')
 
 

figure(3)


plot(Cost3,'r-*')
legend('trace BB') 
 
 

figure(5)
plot(Cost5,'r-*')

 
figure(7)
plot(sum(abs(sum1')),'-*')

figure(8)
plot(sum2,'-*')% 每一次

figure(9)
plot(sum(abs(sum1'))/max(sum(abs(sum1'))),'r-*','LineWidth',1.5)% everage of B guiyihua


figure(10)
plot(number_leader/(iterations),'-*','LineWidth',1.5) % probability
% %set(gca, 'YLim',[0.975 1]); 
% %axis=([0 100  0.975 1]);
% xlabel('Nodes');
% ylabel('Probabitity of selecting as key node');
% %set(gca, 'LineWidth', 1.5);
% set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
% 
% export_fig prob_100_41_dia.eps -painters -transparent
% % 

eig_SS=zeros(N,M);
for i=1:M
 SS=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    SS=SS+expm(-A*(ot*k))*B(:,i)*B(:,i)'*expm(-A'*(ot*k))*ot;
end

eig_SS(:,i)=eig(SS)

end

[N,M]=size(B);
Bb = zeros(N,M);
sum_B=sum(sum1,2);
[max_v order]=sort(sum_B,'descend');%max_v为值，order为在sum_B中的序号

 for i=1:M order(i)=4*i-3; end
for k=1:M
    Bb(order(k),k)=1;
end

WB0=zeros(N,N);

%先求出文章中的WB，即为这里的WB0
for k1=1:tf/ot
    WB0=WB0+expm(-A*(ot*k1))*Bb*Bb'*expm(-A'*(ot*k1))*ot;
end
C=pinv(WB0);
PGME=trace(C)

N-rank(ctrb(A,Bb)) 



% 
%   G=A';
%   G=full(G);
%   num_Node=size(A)
%   figure(13)
%   if num_Node<500 
%    h = view(biograph(G,[],'ShowWeights','off', 'LayoutType', 'equilibrium'));
%    h1 = view(biograph1(G,[],'ShowWeights','on', 'LayoutType', 'equilibrium'));
%    hold on
%    set(h.Nodes,'Shape','circle');
%    set(h.Nodes(Index),'Color',[1 0 0]);
%  
%   set(h.Nodes(Index_final),'Color',[1 0  0 ]); 
%   
%   end
% %   
% 




