close all 
t = tiledlayout(3,4,'TileSpacing','Compact');
load("plotdataMC15t.mat");
% y12all = [];
% y13all = [];
% y23all = [];
% for i = 1:1000
% y12all =  [y12m{i};y12all];
% y13all =  [y13m{i};y13all];
% y23all = [y23m{i};y23all];
% end
y12max = prctile(y12m,90);
y12min = prctile(y12m,10);
y12med = median(y12m);

y12maxX = prctile(y12m,100);
t12 = find(y12maxX<=1,1);
y13max = prctile(y13m,90);
y13min = prctile(y13m,10);
y13med = median(y13m);

y13maxX = prctile(y13m,100);

t13 = find(y13maxX<=1,1);

y23max = prctile(y23m,90);
y23min = prctile(y23m,10);
y23med = median(y23m);

y23maxX = prctile(y23m,100);

t23 = find(y23maxX<1,1);
tm = max(t12,max(t13,t23));
tm*17/60
tspan = (tspan*10)*17/60; 
% figure 
% for i = 1:1000
%     subplot(3,1,1)
%     plot(tspan,y12m(i,:));
%     hold on 
%     subplot(3,1,2)
%     plot(tspan,y13m(i,:));
%     hold on 
%     subplot(3,1,3)
%     plot(tspan,y23m(i,:));
%     hold on 
% end
nexttile(t,[1 3])
for i =1:1000
    L=plot(tspan,y12m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on
f1 =fill([tspan,fliplr(tspan)],[y12min,fliplr(y12max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y12med ,'LineStyle','-','LineWidth',2,'Color','#035AA6');
ylabel('$y_{12}~[\mathrm{\mu rad}]$','Interpreter','latex');
title('Misalignments','Interpreter','latex')
grid on 
hold on
xc = xline(tm*17/60,'LineWidth',1,'Color','k');
%text(tm*17/60,0,'$t^*$','Interpreter','latex');
xticks([0 20 40 tm*17/60 60 80 ])
xticklabels({'0','20','40','$t^*$','60','80'})
set(gca,'TickLabelInterpreter','latex')
axis([0 350*17/60 0 9]);
leg = legend([L p1 xc f1],'Realizations','Median','Worst misalignment $<1 \mu \mathrm{rad}$','10-90 percentile range');
leg.Interpreter = 'latex';
leg.Box = 'off';

nexttile(t)
h1 = histogram(y12m(:,tm),'FaceColor',[126,176,213]/255,'FaceAlpha',0.5);
xlabel('$y_{12}~[\mathrm{\mu rad}]$','Interpreter','latex')
title('Histograms for the measurement at $t^*$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
grid on 
axis([0 1 0 360])
nexttile(t,[1 3])
for i =1:1000
    plot(tspan,y13m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on
f1 =fill([tspan,fliplr(tspan)],[y13min,fliplr(y13max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y13med ,'LineStyle','-','LineWidth',2,'Color','#035AA6');
ylabel('$y_{13}~[\mathrm{\mu rad}]$','Interpreter','latex');
axis([0 350*17/60 0 9]);
grid on 
xc = xline(tm*17/60,'LineWidth',1,'Color','k');
xticks([0 20 40 tm*17/60 60 80 ])
xticklabels({'0','20','40','$t^*$','60','80'})
set(gca,'TickLabelInterpreter','latex')

nexttile(t)
h2 = histogram(y13m(:,tm),'FaceColor',[126,176,213]/255,'FaceAlpha',0.5);
xlabel('$y_{13}~[\mathrm{\mu rad}]$','Interpreter','latex');
grid on 
set(gca,'TickLabelInterpreter','latex')
axis([0 1 0 360])

nexttile(t,[1 3] )
for i =1:1000
    plot(tspan,y23m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on 

f1 =fill([tspan,fliplr(tspan)],[y23min,fliplr(y23max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y23med ,'LineStyle','-','LineWidth',2,'Color','#035AA6');
ylabel('$y_{23}~[\mathrm{\mu rad}]$','Interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
axis([0 350*17/60 0 9]);
xc = xline(tm*17/60,'LineWidth',1,'Color','k');
xticks([0 20 40 tm*17/60 60 80 ])
xticklabels({'0','20','40','$t^*$','60','80'})
grid on
set(gca,'TickLabelInterpreter','latex')

nexttile(t)
h3 = histogram(y23m(:,tm),'FaceColor',[126,176,213]/255,'FaceAlpha',0.5);
xlabel('$y_{23}~[\mathrm{\mu rad}]$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
axis([0 1 0 360])

grid on
% NN1 = h1.NumBins;
% NN2 = h2.NumBins;
% NN3 = h3.NumBins;
% NN = max(NN1,max(NN2,NN3));
% WW1 = h1.BinWidth;
% WW2 = h2.BinWidth;
% WW3 = h3.BinWidth;
% WW = max(NN1,max(NN2,NN3));
% h1.NumBins = NN; 
% h2.NumBins = NN; 
% h3.NumBins = NN; 
% h1.BinWidth = WW;
% h2.BinWidth = WW;
% h3.BinWidth = WW;
% h1.NumBins = NN;
% h2.NumBins = NN;
% h3.NumBins = NN;
% figure 
% histogram(y12m(:,end))
% hold on 
% histogram(y13m(:,end))
% hold on 
% histogram(y23m(:,end))

% load("plotdataMC1000nomomentum.mat");
% y12all = [];
% y13all = [];
% y23all = [];
% for i = 1:MC
% y12all =  [y12m;y12all];
% y13all =  [y13m;y13all];
% y23all = [y23m;y23all];
% end
% y12max = max(y12all);
% y12min = min(y12all);
% y12med = median(y12all);
% 
% y13max = max(y13all);
% y13min = min(y13all);
% y13med = median(y13all);
% 
% y23max = max(y23all);
% y23min = min(y23all);
% y23med = median(y23all);
% 
% nexttile
% f1 =fill([tspan,fliplr(tspan)],[y12min,fliplr(y12max)],[255, 128, 7]/255,'FaceAlpha',0.3);
% hold on 
% p1 = plot(tspan,y12med ,'LineStyle','-','LineWidth',2,'Color',[255, 128, 7]/255);
% 
% nexttile
% f1 =fill([tspan,fliplr(tspan)],[y13min,fliplr(y13max)],[255, 128, 7]/255,'FaceAlpha',0.3);
% hold on 
% p1 = plot(tspan,y13med ,'LineStyle','-','LineWidth',2,'Color',[255, 128, 7]/255);
% 
% nexttile
% f1 =fill([tspan,fliplr(tspan)],[y23min,fliplr(y23max)],[255, 128, 7]/255,'FaceAlpha',0.3);
% hold on 
% p1 = plot(tspan,y23med ,'LineStyle','-','LineWidth',2,'Color',[255, 128, 7]/255);