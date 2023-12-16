clear 
close all
clc

hemo = 1;
load('/data/eggebrecht/data1/Weihao/NeuroDOT_WF/DT/0/metrics_with_0Hz.mat')
LE_dt = LE(:,:,:,hemo).*ROI;
FW_dt = FW(:,:,:,hemo).*ROI;
ER_dt = ER(:,:,:,hemo).*ROI;
LE0_dt = LE_dt(LE_dt>0);
FW0_dt = FW_dt(LE_dt>0 & LE_dt<=8);
ER0_dt = ER_dt(LE_dt>0 & LE_dt<=8);
depth1_dt = depth(roi>0);      % depth for LE and SR
depth2_dt = depth1_dt(LE0_dt<=8);

load('/data/eggebrecht/data1/Weihao/NeuroDOT_WF/MC/0/metrics_with_0Hz.mat')
LE_mc = LE(:,:,:,hemo).*ROI;
FW_mc = FW(:,:,:,hemo).*ROI;
ER_mc = ER(:,:,:,hemo).*ROI;
LE0_mc = LE_mc(LE_mc>0);
FW0_mc = FW_mc(LE_mc>0 & LE_mc<=8);
ER0_mc = ER_mc(LE_mc>0 & LE_mc<=8);
depth1_mc = depth(roi>0);      % depth for LE and SR
depth2_mc = depth1_mc(LE0_mc<=8);    % depth f

% localization error
medianLE_dt = zeros(35,1);
bLE_dt = zeros(35,1);
tLE_dt = zeros(35,1);
medianLE_mc = zeros(35,1);
bLE_mc = zeros(35,1);
tLE_mc = zeros(35,1);
for i=1:1:35
    medianLE_dt(i,1) = median(LE0_dt(depth1_dt(:)<i & depth1_dt(:)>=(i-1)));
    bLE_dt(i,1) = prctile(LE0_dt(depth1_dt(:)<i & depth1_dt(:)>=(i-1)),25,'all');
    tLE_dt(i,1) = prctile(LE0_dt(depth1_dt(:)<i & depth1_dt(:)>=(i-1)),75,'all');
    medianLE_mc(i,1) = median(LE0_mc(depth1_mc(:)<i & depth1_mc(:)>=(i-1)));
    bLE_mc(i,1) = prctile(LE0_mc(depth1_mc(:)<i & depth1_mc(:)>=(i-1)),25,'all');
    tLE_mc(i,1) = prctile(LE0_mc(depth1_mc(:)<i & depth1_mc(:)>=(i-1)),75,'all');
end
fig1 = figure;
plot(0.5:1:34.5,medianLE_dt,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
hold on
px = [0.5:1:34.5,34.5:-1:0.5];
py_dt = [bLE_dt(:,1)',fliplr(tLE_dt(:,1)')];
patch(px,py_dt,[0.9,0.2,0.2],'EdgeColor','none');
alpha(.3);
plot(0.5:1:34.5,medianLE_mc,'o-','MarkerSize',6,'MarkerEdgeColor',[0.2,0.2,0.8],...
            'MarkerFaceColor',[0.2,0.2,0.8],'LineWidth',2,'Color',[0.2,0.2,0.8])
py_mc = [bLE_mc(:,1)',fliplr(tLE_mc(:,1)')];
patch(px,py_mc,[0.2,0.2,0.8],'EdgeColor','none');
alpha(.3);
xlim([0 25])
ylim([0 25])
% ylabel('localization error (mm)')
% xlabel('depth (mm)')
% legend('NIRFASTer','MMCLab')
grid on
set(gca,'FontName','Arial','fontsize',20,'LineWidth',2,'gridlinestyle','--')
set(gcf,'position',[800 800 700 600])

% success rate
SR_dt = zeros(70,1);
SR_mc = zeros(70,1);
for i=1:1:70
    SR_dt(i,1) = 100*sum(LE0_dt(depth1_dt(:)<i/2 & depth1_dt(:)>=(i-1)/2)<=8)/sum(LE0_dt(depth1_dt(:)<i/2 & depth1_dt(:)>=(i-1)/2)>0);
    SR_mc(i,1) = 100*sum(LE0_mc(depth1_mc(:)<i/2 & depth1_mc(:)>=(i-1)/2)<=8)/sum(LE0_mc(depth1_mc(:)<i/2 & depth1_mc(:)>=(i-1)/2)>0);
end
fig2 = figure;
plot(0.25:0.5:34.75,smooth(SR_dt),'.-','MarkerSize',8,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
hold on
plot(0.25:0.5:34.75,smooth(SR_mc),'.-','MarkerSize',8,'MarkerEdgeColor',[0.2,0.2,0.8],...
            'MarkerFaceColor',[0.2,0.2,0.8],'LineWidth',2,'Color',[0.2,0.2,0.8])
xlim([0 30])
ylim([0 100])
% ylabel('success rate %')
% xlabel('depth (mm)')
% legend('NIRFASTer','MMCLab')
grid on
set(gca,'FontName','Arial','fontsize',20,'LineWidth',2,'gridlinestyle','--')
set(gcf,'position',[800 800 700 600])

% FWHM
medianFW_dt = zeros(35,1);
bFW_dt = zeros(35,1);
tFW_dt = zeros(35,1);
medianFW_mc = zeros(35,1);
bFW_mc = zeros(35,1);
tFW_mc = zeros(35,1);
for i=1:1:35
    medianFW_dt(i,1) = median(FW0_dt(depth2_dt(:)<i & depth2_dt(:)>=(i-1)));
    bFW_dt(i,1) = prctile(FW0_dt(depth2_dt(:)<i & depth2_dt(:)>=(i-1)),25,'all');
    tFW_dt(i,1) = prctile(FW0_dt(depth2_dt(:)<i & depth2_dt(:)>=(i-1)),75,'all');
    medianFW_mc(i,1) = median(FW0_mc(depth2_mc(:)<i & depth2_mc(:)>=(i-1)));
    bFW_mc(i,1) = prctile(FW0_mc(depth2_mc(:)<i & depth2_mc(:)>=(i-1)),25,'all');
    tFW_mc(i,1) = prctile(FW0_mc(depth2_mc(:)<i & depth2_mc(:)>=(i-1)),75,'all');
end
bFW_dt(isnan(bFW_dt)) = 0;
tFW_dt(isnan(tFW_dt)) = 0;
bFW_mc(isnan(bFW_mc)) = 0;
tFW_mc(isnan(tFW_mc)) = 0;
fig3 = figure;
plot(0.5:1:34.5,medianFW_dt,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
hold on
px = [0.5:1:34.5,34.5:-1:0.5];
py_dt = [bFW_dt(:,1)',fliplr(tFW_dt(:,1)')];
patch(px,py_dt,[0.9,0.2,0.2],'EdgeColor','none');
alpha(.3);
plot(0.5:1:34.5,medianFW_mc,'o-','MarkerSize',6,'MarkerEdgeColor',[0.2,0.2,0.9],...
            'MarkerFaceColor',[0.2,0.2,0.8],'LineWidth',2,'Color',[0.2,0.2,0.8])
py_mc = [bFW_mc(:,1)',fliplr(tFW_mc(:,1)')];
patch(px,py_mc,[0.2,0.2,0.8],'EdgeColor','none');
alpha(.3);
xlim([0 25])
ylim([5 20])
% ylim([0 25])
% ylabel('FWHM (mm)')
% xlabel('depth (mm)')
% legend('NIRFASTer','MMCLab')
grid on
set(gca,'FontName','Arial','fontsize',20,'LineWidth',2,'gridlinestyle','--')
set(gcf,'position',[800 800 700 600])

% ER
medianER_dt = zeros(35,1);
bER_dt = zeros(35,1);
tER_dt = zeros(35,1);
medianER_mc = zeros(35,1);
bER_mc = zeros(35,1);
tER_mc = zeros(35,1);
for i=1:1:35
    medianER_dt(i,1) = median(ER0_dt(depth2_dt(:)<i & depth2_dt(:)>=(i-1)));
    bER_dt(i,1) = prctile(ER0_dt(depth2_dt(:)<i & depth2_dt(:)>=(i-1)),25,'all');
    tER_dt(i,1) = prctile(ER0_dt(depth2_dt(:)<i & depth2_dt(:)>=(i-1)),75,'all');
    medianER_mc(i,1) = median(ER0_mc(depth2_mc(:)<i & depth2_mc(:)>=(i-1)));
    bER_mc(i,1) = prctile(ER0_mc(depth2_mc(:)<i & depth2_mc(:)>=(i-1)),25,'all');
    tER_mc(i,1) = prctile(ER0_mc(depth2_mc(:)<i & depth2_mc(:)>=(i-1)),75,'all');
end
bER_dt(isnan(bER_dt)) = 0;
tER_dt(isnan(tER_dt)) = 0;
bER_mc(isnan(bER_mc)) = 0;
tER_mc(isnan(tER_mc)) = 0;
fig4 = figure;
plot(0.5:1:34.5,medianER_dt,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
hold on
px = [0.5:1:34.5,34.5:-1:0.5];
py_dt = [bER_dt(:,1)',fliplr(tER_dt(:,1)')];
patch(px,py_dt,[0.9,0.2,0.2],'EdgeColor','none');
alpha(.3);
plot(0.5:1:34.5,medianER_mc,'o-','MarkerSize',6,'MarkerEdgeColor',[0.2,0.2,0.8],...
            'MarkerFaceColor',[0.2,0.2,0.8],'LineWidth',2,'Color',[0.2,0.2,0.8])
py_mc = [bER_mc(:,1)',fliplr(tER_mc(:,1)')];
patch(px,py_mc,[0.2,0.2,0.8],'EdgeColor','none');
alpha(.3);
xlim([0 25])
ylim([5 20])
% ylim([0 25])
% ylabel('Effective resolution (mm)')
% xlabel('depth (mm)')
% legend('NIRFASTer','MMCLab')
grid on
set(gca,'FontName','Arial','fontsize',20,'LineWidth',2,'gridlinestyle','--')
set(gcf,'position',[800 800 700 600])