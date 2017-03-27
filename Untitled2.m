
x = 0:0.01:20;
y1 = 200*exp(-0.05*x).*sin(x);
y2 = 0.8*exp(-0.5*x).*sin(10*x);
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
set(AX(1),'XColor','k','YColor','b');
set(AX(2),'XColor','k','YColor','r');
HH1=get(AX(1),'Ylabel');
set(HH1,'String','Left Y-axis');
set(HH1,'color','b');
HH2=get(AX(2),'Ylabel');
set(HH2,'String','Right Y-axis');
set(HH2,'color','r');
set(H1,'LineStyle','-');
set(H1,'color','b');
set(H2,'LineStyle',':');
set(H2,'color','r');
legend([H1,H2],{'y1 = 200*exp(-0.05*x).*sin(x)';'y2 = 0.8*exp(-0.5*x).*sin(10*x)'});
xlabel('Zero to 20 \musec.');
title('Labeling plotyy');