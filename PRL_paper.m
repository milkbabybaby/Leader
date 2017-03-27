clear all
 clc
 close all
 
M=3;
N=9

%  A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
% A(1,N)=1;

% A good example to show the minimization of  the maximum matching paths
% corresponds to the optimal energy control.
 
%A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
  %A=csvread('E:\network\Electronic networks\s208_st.csvadjacencymatrix.csv');
  
%  A=csvread('E:\network\Food webs\Rhode.csvadjacencymatrix.csv');
 
 A = zeros(N,N);
%   A(2,1) = 1;
%   A(4,2) = 1;
%    A(3,4) = 1;
%    A(5,3) = 1;
%   A(6,5) = 1;
%  A0=randperm(N)';
%  for i=1:N-1
%      A(A0(i+1),A0(i))=1;
%  end 
%      
%  
 
 
  for i=1:N-1
     A(i+1,i)=1
 end
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
 BBB=zeros(N,1)
 aver=1
for iiii=1:aver
    close all
N=size(A)
sum1=zeros(N,M); 
tf=1
N=length(A);
iterations=1
Cost11=zeros(M,N)';
for ii=1:iterations

B=rands(M,N)';
MM=99;


B=sqrt(sqrt(1*(MM-0.5)/trace(B'*B*B'*B)))*B; %initiate the norm to MM-1

Xf=expm(A*tf)*expm(A'*tf);

M=min(size((B)));
WB0=zeros(N,N);
F0=zeros(N,M);  


ot=0.02;
 
 
Iterations_number=1500;
Cost=zeros(Iterations_number,1);
Cost3=zeros(Iterations_number,1);
Cost5=zeros(Iterations_number,1);
v=0.01
 
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



F1=(1/ot)*(2*trace(B'*B*B'*B)-2*MM)*4*B*B'*B; % gradient of norm function
F00=F0/norm(F0);  %unit vector of the gradient of energy function 
F11=F1/norm(F1); % unit vector of the gradient of norm function
F11=-sqrt(sqrt(1*(MM+1)/trace(F11'*F11*F11'*F11)))*F11;    %Normailize gradient of norm function
Cost3(k)=trace(B'*B*B'*B);
 
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
B=sqrt(sqrt(1*(MM+1)/trace(B'*B*B'*B)))*B; %initiate the norm to MM+1
end



%F1=(LL*(4*B*B'*B-4*B));
 

%B=B+0.01*1/(Cost2(k)+1)*(F0-F1);

 

if Cost3(k)>MM ;
else
B=B-v*F22/norm(F22);
end


  k  
end

B=sqrt(sqrt(MM+1))*B/norm(B);
sum1=sum1+(abs(B))/iterations
Cost11=Cost11+B/iterations;
end



% figure(1)
% 
% plot(Cost,'r-*')
%  
%  
% 
% figure(3)
% 
% 
% plot(Cost3,'r-*')
% legend('trace BBBB') 
%  
%  
% 
% figure(5)
% plot(Cost5,'r-*')
% 
% figure(6)
% plot((abs(sum1)),'-*')
%  
% figure(7)
% plot(sum(abs(sum1')),'-*')
% 
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
  
  BB=zeros(N,1)
  for j=1:N
 for i=1:M
     
     BB(j,1)=B(j,i).^2+BB(j,1);
 end 
  end
  
   % plot (BB,'b-*')
    BBB=BB+BBB;
end 
BBBB=BBB/aver
 % plot (A0,BBBB,'b-*')
  
  for i=1:N
  %BBBBB(i,1)=BBBB(A0(i))
  BBBBB(i,1)=BBBB(i)
  end 
  
  plot (BBBBB,'b-*')
  
  
  WB=zeros(N,N);
ot=0.01;
for k=1:tf/ot
    WB=WB+expm(-A*(ot*k))*B*B'*expm(-A'*(ot*k))*ot;
end
eig(WB)
trace(WB)

for i=1:M
  B2=B(:,i)
WB1=zeros(N,N);
ot=0.01;
for k=1:tf/ot
    WB1=WB1+expm(-A*(ot*k))*B2*B2'*expm(-A'*(ot*k))*ot;
end
eig_M(:,i)=eig(WB1)
trace_M(i,1)=trace(WB1)
end

% B4=B(:,2)
% WB2=zeros(N,N);
% ot=0.01;
% for k=1:tf/ot
%     WB2=WB2+expm(-A*(ot*k))*B4*B4'*expm(-A'*(ot*k))*ot;
% end
% eig(WB2)
% trace(WB2)


