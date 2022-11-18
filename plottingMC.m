t = tiledlayout(3,2,'TileSpacing','Compact');
load("plotdataMC1000.mat");
% y12all = [];
% y13all = [];
% y23all = [];
% for i = 1:1000
% y12all =  [y12m{i};y12all];
% y13all =  [y13m{i};y13all];
% y23all = [y23m{i};y23all];
% end
y12max = prctile(y12m,75);
y12min = prctile(y12m,25);
y12med = median(y12m);

y13max = prctile(y13m,75);
y13min = prctile(y13m,25);
y13med = median(y13m);

y23max = prctile(y23m,75);
y23min = prctile(y23m,25);
y23med = median(y23m);
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
nexttile(t,1)
for i =1:1000
    plot(tspan,y12m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on
f1 =fill([tspan,fliplr(tspan)],[y12min,fliplr(y12max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y12med ,'LineStyle','-','LineWidth',2,'Color',[0, 0.4470, 0.7410]);

axis([0 200 0 10]);
nexttile(t,3)
for i =1:1000
    plot(tspan,y13m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on
f1 =fill([tspan,fliplr(tspan)],[y13min,fliplr(y13max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y13med ,'LineStyle','-','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
axis([0 200 0 10]);
nexttile(t, 5 )
for i =1:1000
    plot(tspan,y23m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on 

f1 =fill([tspan,fliplr(tspan)],[y23min,fliplr(y23max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y23med ,'LineStyle','-','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
axis([0 200 0 10]);

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


load("plotdataMC1000nomomentum.mat");
y12max = prctile(y12m,75);
y12min = prctile(y12m,25);
y12med = median(y12m);

y13max = prctile(y13m,75);
y13min = prctile(y13m,25);
y13med = median(y13m);

y23max = prctile(y23m,75);
y23min = prctile(y23m,25);
y23med = median(y23m);

nexttile(t,2)
for i =1:1000
    plot(tspan,y12m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on
f1 =fill([tspan,fliplr(tspan)],[y12min,fliplr(y12max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y12med ,'LineStyle','-','LineWidth',2,'Color',[0, 0.4470, 0.7410]);

axis([0 200 0 10]);
nexttile(t,4)
for i =1:1000
    plot(tspan,y13m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on
f1 =fill([tspan,fliplr(tspan)],[y13min,fliplr(y13max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y13med ,'LineStyle','-','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
axis([0 200 0 155]);
nexttile(t, 6 )
for i =1:1000
    plot(tspan,y23m(i,:),'Color', [.7 .7 .7]);
    hold on
end
hold on 

f1 =fill([tspan,fliplr(tspan)],[y23min,fliplr(y23max)],[126,176,213]/255,'FaceAlpha',0.3);
hold on 
p1 = plot(tspan,y23med ,'LineStyle','-','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
axis([0 200 0 155]);



figure 
thist =  tiledlayout(3,2,'TileSpacing','Compact');
nexttile(thist,2)
histogram(y12m(:,end))
nexttile(thist,4)
histogram(y13m(:,end))
nexttile(thist,6)
histogram(y23m(:,end))
load("plotdataMC1000.mat");
nexttile(thist,1)
histogram(y12m(:,end))
nexttile(thist,3)
histogram(y13m(:,end))
nexttile(thist,5)
histogram(y23m(:,end))


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