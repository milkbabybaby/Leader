clear all;
close all;
clc;


%{
%%%%%%%%%%%%%%%real-life����ֱ������Edge List%%%%%%%%%%%%%
dataset_name = 'Maspalomas'
load_path = strcat('./dataset/',dataset_name,'.mat');
load(load_path);
L = size(A,1)%LΪ�ߵĸ���
N2 = length(unique(A))
[A,~] = Net_Consecutive(double(A)); %��������Ĵ����ֹA��Ų�������������AΪ�ڽӾ���
%}


%%%%%%%%%%%%%%%ER-SF����ֱ�������ڽӾ���%%%%%%%%%%%%%%%%%%%
%   A=ErSfNetGen(2,100,5);
%   dataset_name = 'A100ER'
%   L=sum(sum(A))

%load('D:\code\dataset\Foodweb\Rhode.mat')
load('D:\code\Leader\A.mat');

% L = size(A,1)%LΪ�ߵĸ���
% N2 = length(unique(A))
% [A,~] = Net_Consecutive(double(A)); 
N=size(A,1);
A=A';

%A=[ 0 1 0; 0 0 1; 0 0 0];
%{
%����stem���ڽӾ���
A=[0 0 0 0 0 1 0 1 0 0;
   0 0 0 0 0 0 0 0 0 0;
   0 1 0 0 1 1 0 1 0 0;
   1 0 0 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 1 0;
   0 1 0 0 0 0 0 1 0 0;
   0 0 0 0 0 0 1 0 0 0;
   0 0 0 0 1 1 0 1 0 0;
   0 0 0 1 0 0 0 1 1 0;];

%}
N=size(A,1)
M=60;

sum1=zeros(N,M);
v=0.02; %��������
ot=0.02; %����ʱ�䲽
iterations=1;  %�ܵ�ѭ������
Iterations_number=500; 
tf = 1;


BB=zeros(N,M);

for ii=1:iterations
    
    B=randn(N,M); %BΪ����Ȩ�ؾ��������ʼֵ
    MM=M; %�涨��һ��normֵ
    
    B=(sqrt(1*(MM+0.1)/trace(B'*B)))*B; %initiate the norm to MM-1

    M=min(size((B)));  %M����һ������¾���M
    WB0=zeros(N,N); 
    F0=zeros(N,M);
    Xf=expm(A*tf)*expm(A'*tf); %�Ծ�����ָ������ÿ��Ԫ����ָ��
    
    Cost=zeros(Iterations_number,1);
    Cost3=zeros(Iterations_number,1);
    Cost5=zeros(Iterations_number,1);

    
    %1200���ݶ��½�����
    for k=1:Iterations_number
        fprintf('No. %d Epoch, No. %d Descent \n',ii,k);
        
        WB0=zeros(N,N);
        F0=zeros(N,M);
        %����������е�WB����Ϊ�����WB0
        for k1=1:tf/ot
            WB0 = WB0+expm(A*(ot*k1))*B*B'*expm(A'*(ot*k1))*ot;
        end
        C = pinv(WB0);
        
        Cost(k)=trace(C*Xf);
        Cost3(k)=trace(B'*B);
        
        %����������е�����������B������ݶȡ�
        for k1=1:tf/ot
            F0=F0-(expm(A'*(ot*k1))*C*Xf*C*expm(A*(ot*k1)))*B*ot;
        end
        F00=F0/norm(F0); 
        F22=F0; %F22Ϊ���ε������ݶ�
        
        %F1ΪNorm������B���ݶ�
        F1=(1/ot)*(2*trace(B'*B)-2*MM)*2*B; % gradient of norm function
        F11=F1/norm(F1); %��һ��Norm�����ݶ�
        F11=-(sqrt(1*(MM+0.001)/trace(F11'*F11)))*F11;    %Normailize gradient of norm function
        
        Cost5(k)=-sum(F00(:).*F11(:))/(norm(F00(:))*norm(F11(:)));
        
        %���cost5�ӽ�-1�������ݶ��½�����ѭ��
        if  Cost5(k)+1<0.000005 
            break;
        end
        
        if  Cost3(k)>=MM  %MM is the fixed energy bound
            %F0=reshape((eye(N*M)-F1(:)*pinv(F1(:)))*F0(:),N,M);  %the direction of gradient energy function  (projection)
            
            %��һ������ȫΪ�˽�ʡ�ռ��������������
            F = F1(:); %FΪ������NM*1
            P = F'*F;
            P = 1./P;
            Fp = P*F'; %FpΪ������1*NM
            
            FF = F0(:); %FFΪ������NM*1
            Fv = zeros(N*M,1);
            C = Fp*FF;
            Fv = FF - F.*C;
            F0 = reshape(Fv,N,M);%the direction of gradient energy function  (projection)
     
            F0=F0/norm(F0); %unit vector of gradient energy function         
            if sum(sum(F0.*F22))<=0
                F0=-F0;
            end
  
            B = B-v*F0/norm(F0);  % F0 is the projected direction
            B = sqrt(1*(MM+0.001)/trace(B'*B))*B; %initiate the norm to MM+1
        else
            B = B-v*F0/norm(F0);
        end
   
    end
    
    BB=BB+abs(B);
end

%���������logͼ
figure(1)
plot(log10(Cost),'r-*')

% figure(2)
% plot(Cost3,'r-*')
% legend('trace BB')

figure(3)
plot(Cost5,'r-*')


Init = Cost(1)
PGM = Cost(end)


%������PGME
%����Science paper�ķ�����Bs
[N,M]=size(B);
Bb = zeros(N,M);
sum_B=sum(abs(BB),2);
[max_v order]=sort(sum_B,'descend');%max_vΪֵ��orderΪ��sum_B�е����
for k=1:M
    Bb(order(k),k)=1;
end

WB0=zeros(N,N);

%����������е�WB����Ϊ�����WB0
for k1=1:tf/ot
    WB0=WB0+expm(A*(ot*k1))*Bb*Bb'*expm(A'*(ot*k1))*ot;
end
C=inv(WB0);
PGME=trace(C*Xf)

N-rank(ctrb(A,Bb)) % ��zhiΪ0

 save('D:\code\A_conference','Bb')

%save_path = strcat('PGM-M',num2str(M),'tf1_',dataset_name);
%save(save_path,'Cost','B','Bb','Init','PGM', 'PGME','Iterations_number');
%save(save_path,'Cost','B','Bs1','Bs2','Init','PGM', 'PGME1','PGME2','Iterations_number');

