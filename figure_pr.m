clear;
clc;
x=[ 10, 100, 1000]
y1=[0.7767,0.7930,0.7960]
y2=[0.9171,0.9867,0.9985]
y3=[0.9035,0.9236,0.9242]
y4=[0.7636,0.7689,0.77]

semilogx(x,y1,'g-*','LineWidth',1.5,'markersize',8) % probability
hold on
semilogx(x,y2,'b-.x','LineWidth',1.5,'markersize',8) % probability
hold on
semilogx(x,y3,'r:^','LineWidth',1.5,'markersize',8) % probability
hold on
semilogx(x,y4,'m--o','LineWidth',1.5,'markersize',8) % probability
%set(gca, 'YLim',[0.975 1]); 
%axis=([0 100  0.975 1]);
xlabel('Number of experiments n','FontName','Arial','FontSize',10,'FontWeight','Bold');
ylabel('The  correlation  coefficient \rho','FontName','Arial','FontSize',10,'FontWeight','Bold');
%set(gca, 'LineWidth', 1.5);
legend('Network I','Network II','Network III ','Network IV ','FontName','Arial','FontSize',10,'FontWeight','Bold','Location', 'East' ) 

set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')

export_fig pr.eps -painters -transparent