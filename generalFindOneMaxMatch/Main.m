%Main������Ѱ��һ��MMƥ��
%����EdgeList��Ȼ�󷵻�DriverNode���
%A_12��1��2
clear all;
clc;
%digits(80)
%N=8

%A=BAgraph_dir_SF(100,15,6)
A=rand_ER(100,0.15)
  %A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
  %A=A'

%%%%%%%%%%%%%%real-life����ֱ������Edge List%%%%%%%%%%%%%
%dataset_name = 'celegans';
%load_path = strcat('./dataset/',dataset_name,'.mat');

%load('D:\code\dataset\circuit-s208.mat');
%load('D:\code\dataset\CrystalC.mat')
%Rhode, Michigan, cons-frequency-rev    cons-quality-rev    CrystalC   CrystalID  gramwet  Mondego  
%  Florida   StMarks
% load('A')
%  A=A';
%[A,~] = Net_Consecutive(double(A)); %��������Ĵ����ֹA��Ų�������������AΪ�ڽӾ���,��Ϊ�����Ϊ�յ�  A_12 shi 1 dao 2
%

% %%%%%%%%%%%%%%%ER-SF����ֱ�������ڽӾ���%%%%%%%%%%%%%%%%%%%
 %A=ErSfNetGen(2,200,1);%kkk=1ΪER���磬����ΪSF��BA������
% dataset_name = 'BA200-u8';
% L=sum(sum(A))
% save_path = strcat(dataset_name);




G = generalNewBG(A);
MM_mat = generalFindOneMaxMatch(G);%����ƥ�����
if issparse(MM_mat)%����ϡ�����
    MM_mat = full(MM_mat);%ת��Ϊ��ͨ����
end
N = size(MM_mat);%�ܵĽڵ���
Parent_num = sum(MM_mat,2);%�ҳ�ÿ���ڵ�ĸ��ڵ�����ƥ����������
Driver_Nodes = [];
for i=1:N
    if Parent_num(i)==0%û�и��ڵ�ľ���Driver_Nodes
        Driver_Nodes = [Driver_Nodes i];
    end
end
% 
% save_path = strcat(dataset_name,'MM')
% save(save_path,'MM_mat','Driver_Nodes')
 save('A','A','Driver_Nodes');
 
 
N=size(A,1)
ND=size(Driver_Nodes,2);
Bb = zeros(N,ND);
for k=1:ND
    Bb(Driver_Nodes(k),k)=1;
end
N-rank(ctrb(A',Bb),0)




for i=1:length(Driver_Nodes)

p=Driver_Nodes(i);
path=[p]
row_nodes=MM_mat(:,p);  %match�ߵ����Ϊp,p��j������
 while length(find(row_nodes==1))
     p=find(row_nodes==1)
     path=[path,p]
row_nodes=MM_mat(:,p);
 end
 path_driver(i,1)={path};
 end 
u=2*sum(sum(A))/N     
     diag(A);
% kin=sum(A,1)
% figure
%  Nbin1 = histc(kin,min(kin):max(kin)+1);
%     Nbin1=Nbin1(1:length(Nbin1)-1)/N
% subplot(1,2,1)
% bar(min(kin):max(kin),Nbin1);xlim([-1, max(kin)+1]);
% text(0.85,0.92,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontName','Arial','FontSize',10,'FontWeight','Bold')
% xlabel('degree of all nodes','FontName','Arial','FontSize',10,'FontWeight','Bold');
% ylabel('Probability','FontName','Arial','FontSize',10,'FontWeight','Bold');
% %set(gca, 'LineWidth', 1.5);
% set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
