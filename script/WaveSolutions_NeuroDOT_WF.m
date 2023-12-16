clc
clear
close all


%% Compare TPSF
frequency = 0;     % Hz
path = '/data/eggebrecht/data1/Weihao/NeuroDOT_WF';
load([path,'/mesh.mat'])
load([path,'/MC/flux_wavelength2_source20.mat'])

m = 138;
mc_tpsf = flux.data(mesh.elements(mesh.source.int_func(m,1),1),:);
mc_time = 0:1e-10:5e-9-1e-10;
src_coor = mesh.source.coord(20,:);
det_coor = mesh.nodes(mesh.elements(mesh.source.int_func(m,1),1),:);

% plot pulse signal
figure
semilogy(mc_time,abs(mc_tpsf(1,:)),'o-','LineWidth',1.5)
% xlim([0 4e-9])
% ylim([1e-13 1e4])
ti = title(['Rsd = ',num2str(floor(pdist2(src_coor,det_coor))),'mm'],'Interpreter','latex');
ti.FontSize = 19;
xl = xlabel('time /s','Interpreter','latex');
xl.FontSize = 17;
yl = ylabel('Log10 TPSF /mm$^{-2}$s$^{-1}$','Interpreter','latex');
yl.FontSize = 17;

