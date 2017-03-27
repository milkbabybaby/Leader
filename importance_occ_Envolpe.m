clc;
clear;

%% figure plot stem
load('50_10_1000_stem.mat')

figure(1)
subplot(2,1,1)
b0=sum(abs(sum1'))/max(sum(abs(sum1')));
t=1:length(b0);
plot(sum(abs(sum1'))/max(sum(abs(sum1'))),'r-.*','LineWidth',1.5);% everage of B guiyihua
%set(gca, 'YLim',[0.9 1.08]); 
set(gca, 'YLim',[0 1.1]); 
%axis=([0 100  0.5  1]);
xlabel('Node index');
ylabel('r');
%set(gca, 'LineWidth', 1.5);
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
% [up,lo] = envelope(b0,15,'peak');
[up,lo] = envelope(b0,15,'peak');
hold on
plot(t,up,'b:',t,lo,'b:','linewidth',1.5)
hold off


% grid on 
% grid minor
%export_fig importance_100_cycle.eps -painters -transparent
subplot(2,1,2)
b0=number_leader/(iterations);
t=1:length(b0);
plot(number_leader/(iterations),'b-.*','LineWidth',1.5) % probability
set(gca, 'YLim',[-0.1 1.1]); 
%axis=([0 50 0 1]);
xlabel('Node index');
ylabel('p');
b0(6)=0.9
%set(gca, 'LineWidth', 1.5);
[up,lo] = envelope(b0,7,'peak');
hold on
plot(t,up,'r:',t,lo,'r:','linewidth',1.5)
hold off
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')

samexaxis('xmt','on','ytac','join','yld',1)
%% figure plot cycle

load('100_25_1000_CYCLE.mat')
figure(2)
subplot(2,1,1)
plot(sum(abs(sum1'))/max(sum(abs(sum1'))),'r--*','LineWidth',1.5);% everage of B guiyihua
set(gca, 'YLim',[0.9 1.08]); 
%axis=([0 100  0.5  1]);
xlabel('Node index');
ylabel('r');
%set(gca, 'LineWidth', 1.5);
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
% grid on 
% grid minor
%export_fig importance_100_cycle.eps -painters -transparent
% 
subplot(2,1,2)
plot(number_leader/(iterations),'b--*','LineWidth',1.5) % probability
set(gca, 'YLim',[0.0 0.55]); 
%axis=([0 100  0.0 0.5]);
xlabel('Node index');
ylabel('p');
%set(gca, 'LineWidth', 1.5);
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
samexaxis('xmt','on','ytac','join','yld',1)
%export_fig prob_100_cycle.eps -painters -transparent


%% plot cycle
% load('100_25_1000_CYCLE.mat')
% 
% figure(1)
% plot(sum(abs(sum1'))/max(sum(abs(sum1'))),'r-*','LineWidth',1.5);% everage of B guiyihua
% set(gca, 'YLim',[0.9 1.08]); 
% %axis=([0 100  0.5  1]);
% xlabel('Node index');
% ylabel('r');
% %set(gca, 'LineWidth', 1.5);
% set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
% % grid on 
% % grid minor
% export_fig importance_100_cycle.eps -painters -transparent
% % 
% figure(2)
% plot(number_leader/(iterations),'-*','LineWidth',1.5) % probability
% set(gca, 'YLim',[0.0 0.5]); 
% %axis=([0 100  0.0 0.5]);
% xlabel('Node index');
% ylabel('p');
% %set(gca, 'LineWidth', 1.5);
% set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')



%% plot dialation
load('100_25_41_1000_dia.mat')
figure(3)
subplot(2,1,1)
b0=sum(abs(sum1'))/max(sum(abs(sum1')));
t=1:length(b0);
plot(sum(abs(sum1'))/max(sum(abs(sum1'))),'r-.*','LineWidth',1.5);% everage of B guiyihua
set(gca, 'YLim',[0.1 1]); 
%set(gca, 'YLim',[0.9 1.08]); 
%axis=([0 100  0.5  1]);
xlabel('Node index');
ylabel('r');
%set(gca, 'LineWidth', 1.5);
[up,lo] = envelope(b0,4,'peak');
hold on
plot(t,up,'b:',t,lo,'b:','linewidth',1.5)
hold off
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
% grid on 
% grid minor
%export_fig importance_100_cycle.eps -painters -transparent
% 
subplot(2,1,2)
b0=number_leader/(iterations);
t=1:length(b0);
plot(number_leader/(iterations),'b-.*','LineWidth',1.5) % probability
set(gca, 'YLim',[-0.1 1.1]); 
%set(gca, 'YLim',[0  1.1]); 
axis=([0 100  0.0 1]);
xlabel('Node index');
ylabel('p');
%linkaxes(h,'x');
[up,lo] = envelope(b0,3,'peak');
hold on
plot(t,up,'r:',t,lo,'r:','linewidth',1.5)
hold off
set(gca, 'LineWidth', 1.5);
set(gca,'LineWidth', 1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
samexaxis('xmt','on','ytac','join','yld',1)
%export_fig prob_100_cycle0.eps -painters -transparent
