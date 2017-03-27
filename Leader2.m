clear all
 clc
 close all
 

N=100
M=25
NN=41
%  A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
% A(1,N)=1;

% A good example to show the minimization of  the maximum matching paths
% corresponds to the optimal energy control.
 
%A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
  %A=csvread('E:\network\Electronic networks\s208_st.csvadjacencymatrix.csv');
  
%  A=csvread('E:\network\Food webs\Rhode.csvadjacencymatrix.csv');
 
 A = zeros(N,N);
%  A(2,1)=1
% A(3,2)=1
% A(4,3)=1
% A(5,4)=1
% A(6,5)=1
% A(7,3)=1
% A(8,7)=1
% A(9,8)=1
% A(10,9)=1
% A(7,10)=1
% A(11,1)=1
% A(12,11)=1
% A(13,11)=1
% A(14,13)=1
% 
% 
% A(15,4)=1
% A(16,15)=1
% A(17,16)=1
% A(18,17)=1
% A(15,18)=1






%  A(2,1) = 1;
%  A(3,2) = 1;
%   A(4,3) = 1;
%   A(5,4) = 1;
%  A(6,5) = 1;

  for i=1:N-1
     A(i+1,i)=1
  end
%  
% % A(1,N)=1 %cycle
% 
%  A(NN,1)=1 %dialtion
%  A(NN,NN-1)=0
%  

% A(5,4) = 1;
% A(6,5) = 1;
% A(7,6) = 1;
% A(8,7) = 1;
% A(9,3) = 1;
% A(10,9) = 1;
% A(11,10) = 1;

%A(N,1)=1;

% A=[0 1 0 0 0 0 0 0  0 0 0 0;
%    0 0 1 0 0 0 1 0  0 0 0  0;
%    0 0 0 1 0 0 0 0  0 0 0  0;
%    0 0 0 0 1 0 0 0  0 0 0  0;
%    0 0 0 0 0 1 0 0  0 0 0  0;
%    0 0 0 0 0 0 0 0  0 0 0  0;
%    0 0 0 0 0 0 0 1  0 0 0  0;
%    0 0 0 0 0 0 0 0  1 0 0  0;
%    0 0 0 0 0 0 0 0  0 1 0  0;
%    0 0 0 0 0 0 0 0  0 0 1  0;
%    0 0 0 0 0 0 0 0  0 0 0  1
%    0 0 0 0 0 0 0 0  0 0 0 0];
% 
% 
% 
%  A=A';

N=size(A)
sum1=zeros(N,M); 
sum2=zeros(N,1); 
sum3=zeros(N,1);
tf=2
N=length(A);

Cost11=zeros(M,N)';
iterations=100
for ii=1:iterations

B=rands(M,N)';
MM=99;


B=(sqrt(1*(MM-0.5)/trace(B'*B)))*B; %initiate the norm to MM-1

Xf=expm(A*tf)*expm(A'*tf);

M=min(size((B)));
WB0=zeros(N,N);
F0=zeros(N,M);  


ot=0.0025;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=inv(WB0); 

for k=1:tf/ot
 F0=F0+expm(A'*(ot*k))*C*Xf*C*expm(A*(ot*k))*B*ot;
end
 
 
Iterations_number=1500;
Cost=zeros(Iterations_number,1);
Cost3=zeros(Iterations_number,1);
Cost5=zeros(Iterations_number,1);
v=0.006
 
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

  F22=F0;   
 
Cost(k)=trace(C*Xf);


%F1=(1/ot)*((2*trace(B'*B*B'*B)*(4*B*B'*B-4*B)-2*MM*(4*B*B'*B-4*B)));



F1=(1/ot)*(2*trace(B'*B)-2*MM)*2*B; % gradient of norm function
F00=F0/norm(F0);  %unit vector of the gradient of energy function 
F11=F1/norm(F1); % unit vector of the gradient of norm function
F11=-(sqrt(1*(MM+1)/trace(F11'*F11)))*F11;    %Normailize gradient of norm function
Cost3(k)=trace(B'*B);
 
Cost5(k)=-sum(F00(:).*F11(:))/(norm(F00(:))*norm(F11(:)));

if  Cost5(k)+1<0.0005
    break
end


if  Cost3(k)>=MM  %MM is the fixed energy bound 
F0=reshape((eye(N*M)-F1(:)*pinv(F1(:)))*F0(:),N,M);  %the direction of gradient energy function  (projection) 
F0=F0/norm(F0); %unit vector of gradient energy function 

if sum(sum(F0.*F22))<=0
    F0=-F0;
end



B=B-v*F0/norm(F0);  % F0 is the projected direction 
B=(sqrt(1*(MM+1)/trace(B'*B)))*B; %initiate the norm to MM+1
end



%F1=(LL*(4*B*B'*B-4*B));
 

%B=B+0.01*1/(Cost2(k)+1)*(F0-F1);

 

if Cost3(k)>MM ;
else
B=B-v*F22/norm(F22);
end

k
   
end
BBB_tensor(:,:,ii)=B;
 save('lalalala')
 
 
ii
B=(sqrt(MM+1))*B/norm(B);
sum1=sum1+(abs(B))/iterations;                                   
BB=sum((abs(B)),2);

Cost11=Cost11+B/iterations;

[YY,pos]=sort(BB,'descend');
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

figure(6)
plot((abs(sum1)),'-*')
 
figure(7)
plot(sum(abs(sum1')),'-*')


figure(9)
plot(sum(abs(sum1'))/max(sum(abs(sum1'))),'r-*','LineWidth',1.5)% everage of B guiyihua
%set(gca, 'YLim',[0.975 1]); 
%axis=([0 100  0.975 1]);
xlabel('Nodes');
ylabel('importance index');
%set(gca, 'LineWidth', 1.5);
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
% grid on 
% grid minor
export_fig importance_100_41_dia.eps -painters -transparent

figure(10)
plot(number_leader/(iterations),'-*','LineWidth',1.5) % probability
%set(gca, 'YLim',[0.975 1]); 
%axis=([0 100  0.975 1]);
xlabel('Nodes');
ylabel('Occurrence rate of selecting as key node');
%set(gca, 'LineWidth', 1.5);
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')

export_fig prob_100_41_dia.eps -painters -transparent
% 
  G=A';
  G=full(G);
  num_Node=size(A)
  figure(13)
  if num_Node<500 
   h = view(biograph(G,[],'ShowWeights','off', 'LayoutType', 'equilibrium'));
   h1 = view(biograph1(G,[],'ShowWeights','on', 'LayoutType', 'equilibrium'));
   hold on
   set(h.Nodes,'Shape','circle');
   set(h.Nodes(Index),'Color',[1 0 0]);
 
  set(h.Nodes(Index_final),'Color',[1 0  0 ]); 
  
  end
%   
  BB=zeros(N,1)
  for j=1:N
 for i=1:M
     
     BB(j,1)=B(j,i).^2+BB(j,1);
 end 
  end
  prop=number_leader/(iterations);
 idx_notleader=find(prop<0.7);
 prop_notleader= prop(idx_notleader);
 for i=1:M
    
     prob_group(i,1)=sum(prop_notleader((N/M-1)*(i-1)+1:(N/M-1)*i))/(N/M-1);
     idx_group(i,1)=sum(idx_notleader((N/M-1)*(i-1)+1:(N/M-1)*i))/(N/M-1);
    
 end
 figure(11)
 plot(idx_group,prob_group)
 
 lin=[]
 for kkk=1:25
 lin=[lin 1 0 0 0 0]
 end 
 
 im=sum(abs(sum1'))/max(sum(abs(sum1')));
 prob=number_leader/(iterations);
  im*prob/(norm(im)*norm(prob))
  
  var(sum(abs(sum1'))/max(sum(abs(sum1'))))
  var(number_leader/(iterations))
  
  
  
 B=zeros(N,M);
nn=N/M
for i=1:M
    B(nn*i-nn+1,i)=1;
end  

 WB0=zeros(N,N);
for k1=1:tf/ot
    WB0=WB0+expm(A*(ot*k1))*B*B'*expm(A'*(ot*k1))*ot;
end
 C=pinv(WB0); 
 
trace(C*Xf)


 

for jjj =1 :5
    for jj=1:25
        for j=1:100
         
    BBB_log(j,jj,jjj)=exp(abs(BBB_tensor(j,jj,jjj)));
        end
    end
    BBB_logsum(:,jjj)=sum(abs(BBB_log(:,:,jjj)),2);
end
BBB_sum=sum(BBB_logsum,2)
[YY,poss]=sort(BBB_sum,'descend');

syms t
eA=expm(A*t)
eAt=subs(eA,t,2)
[v d]=eig(eAt)



 ind=randperm(N-1,M-1)+1
 


%  
%  
%  B=zeros(8,4);
%  B(1,1)=1;
%  B(2,2)=1;
%  B(5,3)=1;
%  B(1,4)=1;
%  
  B=zeros(N,M);
ind=[2]
 B(1,1)=1;
 for i = 1:M-1
     B(ind(i),i+1)=1;
 end 

 
 
 rank(ctrb(A,B))
 
 syms t
 W0= int(expm(-A*t)*B*B'*expm(-A'*t),t,0,tf)
  
  W0=zeros(N,N);
for k1=1:tf/ot
    W0=W0+expm(-A*(ot*k1))*B*B'*expm(-A'*(ot*k1))*ot;
end
 
WW0=inv(W0);
 trace(inv(W0))
 
 
 [YYY,poss]=sort(sum(abs(sum1')),'descend');
  poss(1:M)
 
  
  sum1=zeros(N,M); 
for i=1:100
       B=BBB_tensor(:,:,i);
  sum1=sum1+(abs(B))/100;    
  
end 


  for kk=1:N
e=find(nodes_leader(:,1:100)==kk); 
number_leader(kk,1)=length(e);
  end 


   
 im=sum(abs(sum1'))/max(sum(abs(sum1')));
 prob=number_leader/(100);
  im*prob/(norm(im)*norm(prob))
 