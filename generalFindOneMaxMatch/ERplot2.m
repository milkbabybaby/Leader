%%% control cost between r-pgme and rsm for ER networks N=100 \deltaM=30
clc
clear
load 'ERplot00.mat'
x=C(1,:)';
y1=C(2,:)';
y2=C(3,:)';
plot(x,y1,'rx','LineWidth',1.5,'markersize',10)
hold on 
plot(x,y2,'b+','LineWidth',1.5,'markersize',8)
p=polyfit(x,y1,5); 
xx=linspace(min(x),max(x));  % 绘图用到的点的横坐标 
yy=polyval(p,xx);   % 拟合曲线的纵坐标 
plot(xx,yy,'c','LineWidth',4); 
hold on 

p=polyfit(x,y2,5); 
xx=linspace(min(x),max(x));  % 绘图用到的点的横坐标 
yy=polyval(p,xx);   % 拟合曲线的纵坐标 
plot(xx,yy,'g','LineWidth',4); 

xlabel('<k>/2','FontName','Arial','FontSize',10,'FontWeight','Bold');
ylabel('lg E(B)','FontName','Arial','FontSize',10,'FontWeight','Bold');
%set(gca, 'LineWidth', 1.5);
legend('R-PGME','RSM','Fitting line for R-PGME','Fitting line for RSM' ) 
%set(gca, 'XLim',[0 32]); 
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')

export_fig er.eps -painters -transparent

