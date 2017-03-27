%%% control cost between r-pgme and rsm for SF networks N=100 M=40

clc
clear
load 'cmean0111.mat'
x=Cmean(1,:)';
y1=Cmean(2,:)';
y2=Cmean(3,:)';
plot(x,y1,'rx','LineWidth',1.5,'markersize',10)
hold on 

plot(x,y2,'b+','LineWidth',1.5,'markersize',10)
hold on 
% p=polyfit(x,y1,5); 
% xx=linspace(min(x),max(x));  % 绘图用到的点的横坐标 
% yy=polyval(p,xx);   % 拟合曲线的纵坐标 
% plot(xx,yy,'c','LineWidth',3); 
% hold on 
p=polyfit(x,y1,6); 
%p=  [  9.177e-06  -0.0008788      0.03236      -0.5655      7.627  ]
xx=linspace(min(x),max(x));  % 绘图用到的点的横坐标 
yy=polyval(p,xx);   % 拟合曲线的纵坐标 
plot(xx,yy,'c','LineWidth',4); 
 hold on 

% p=polyfit(x,y2,2); 
 %p=[0.0003319,-0.01772, 0.3464, -2.958 ,13.92]
  % p=[  2.074e-05  -0.002215       0.0866      -1.479       13.92  ]
     a =       17.39  %(-39.02, 73.8)
       b =     -0.5185  %(-1.322, 0.2851)
       c =       4.427 % (2.381, 6.473)
       d =    0.001449  %(-0.03097, 0.03387)
%     a =       17.39  
%        b =     -0.2593  
%        c =       4.427  
%        d =   0.0007247 
    xx=linspace(min(x),max(x)); 
         yy= a*exp(b*xx) + c*exp(d*xx)

     
 % 绘图用到的点的横坐标 
 %yy=polyval(p,xx);   % 拟合曲线的纵坐标 
 plot(xx,yy,'g','LineWidth',4); 

xlabel('<k>/2','FontName','Arial','FontSize',10,'FontWeight','Bold');
ylabel('lg E(B)','FontName','Arial','FontSize',10,'FontWeight','Bold');
%set(gca, 'LineWidth', 1.5);
legend('R-PGME','RSM','Fitting line for R-PGME','Fitting line for RSM') 

set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
%set(gca, 'XLim',[8 32]); 

export_fig cmean.eps -painters -transparent

% figure
% 
% plot(log10(x),y1)