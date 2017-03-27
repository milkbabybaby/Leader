clc;
clear;
load('100_25_1000_cycle.mat')

x = 1:1:100
y1 = sum(abs(sum1'))/max(sum(abs(sum1')));
y2 = number_leader/iterations;
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
set(AX(1),'XColor','k','YColor','r','Ylim',[0.5,1]);
set(AX(2),'XColor','k','YColor','b','Ylim',[0 0.5]);

HH1=get(AX(1),'Ylabel');
set(HH1,'String','Importance index','LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold');
set(HH1,'color','r');
HH2=get(AX(2),'Ylabel');
set(HH2,'String','Occurrence rate of selecting as key node','LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold');
set(HH2,'color','b');

set(H1,'LineStyle','-','color','r','Marker','*','LineWidth', 1.5);

set(H2,'LineStyle',':','color','b','Marker','*','LineWidth', 1.5);

legend([H1,H2],{'Importance index','Occurrence rate'},0);
%set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
xlabel('Nodes');
%title('Labeling plotyy');
