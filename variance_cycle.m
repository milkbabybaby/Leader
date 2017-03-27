clear;
clc;
close all;
x=[1, 10, 100, 1000]
y1=[0.0089,9.9683*10^-4,1.2078*10^-4,1.5656*10^-5]
y2=[0.1894, 0.0132, 0.0019, 2.1226*10^-4]
   
   
   
figure
   semilogx(x,y1,'r*:','LineWidth', 1.5)
   hold on 
    semilogx(x,y2,'bo-','LineWidth', 1.5)
   
      %set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('experiment times n','FontName','Times New Roman','FontWeight','bold');
 ylabel('variance','FontName','Times New Roman','FontWeight','bold');
 legend('m','p','FontName','Times New Roman','FontWeight','bold')
 set(gca,'LineWidth', 1.5,'FontName','Times New Roman','FontSize',10,'FontWeight','Bold')
  export_fig variance_cycle.eps -painters -transparent 
   
 