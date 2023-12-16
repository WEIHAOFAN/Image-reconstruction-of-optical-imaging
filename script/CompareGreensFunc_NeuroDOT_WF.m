clc
close all
clear 

frequency = 400;     % MHz
src = 20;
path = '/data/eggebrecht/data1/Weihao/NeuroDOT_WF';
load([path,'/mesh.mat'])
load([path,'/MC/',num2str(frequency),'/GFunc_MC_',num2str(frequency),'MHz_.mat'])   % MMCLab
load([path,'/DT/',num2str(frequency),'/GFunc_LUMO_36tiles_on_atlas_testFD.mat'])   % NIRFASTer

%% Process the Green's functions
% Get Green's function by MMCLab
fooRmc = real(squeeze(mcGs(2,src,:)));
fooImc = imag(squeeze(mcGs(2,src,:)));
fooAmc = sqrt(fooRmc.^2+fooImc.^2);
gs1mc = fooAmc;
gs1mc = gs1mc./max(gs1mc);
m0mc.nodes = mesh.nodes;
m0mc.elements = mesh.elements;
m0mc.data = gs1mc;
gsIdxmc = find(gs1mc==max(gs1mc(:)));
m1mc = CutMesh(m0mc,find(m0mc.nodes(:,1)>mesh.nodes(gsIdxmc,1)));
[Iamc,Ibmc] = ismember(m1mc.nodes,m0mc.nodes,'rows');
Ibmc(Ibmc==0) = [];
m1mc.region = mesh.region(Ibmc);
m1mc.data = m0mc.data(Ibmc);
m1mc.data(m1mc.data<0) = 0;
m1mc.data = log10(10^(10)*m1mc.data);
foomc = atan2(fooImc,fooRmc).*(gs1mc>1e-10);
m2mc = m1mc;
m2mc.data = -foomc(Ibmc);
for i=1:length(m2mc.data)
    if m2mc.data(i,1)>0
        m2mc.data(i,1) = m2mc.data(i,1)-2*pi;
    end
end
m1mc.data = -m1mc.data;
m2mc.data = -m2mc.data;

% Get Green's function by NIRFASTer
fooRdt = real(squeeze(Gs(2,src,:)));
fooIdt = imag(squeeze(Gs(2,src,:)));
fooAdt = sqrt(fooRdt.^2+fooIdt.^2);
gs1dt = fooAdt;
gs1dt = gs1dt./max(gs1dt);
m0dt.nodes = mesh.nodes;
m0dt.elements = mesh.elements;
m0dt.data = gs1dt;
gsIdxdt = find(gs1dt==max(gs1dt(:)));
m1dt = CutMesh(m0dt,find(m0dt.nodes(:,1)>mesh.nodes(gsIdxmc,1)));
[Iadt,Ibdt] = ismember(m1dt.nodes,m0dt.nodes,'rows');
Ibdt(Ibdt==0) = [];
m1dt.region = mesh.region(Ibdt);
m1dt.data = m0dt.data(Ibdt);
m1dt.data(m1dt.data<0) = 0;
m1dt.data = log10(10^(10)*m1dt.data);
foodt = atan2(fooIdt,fooRdt).*(gs1dt>1e-10);
m2dt = m1dt;
m2dt.data = foodt(Ibdt);
for i=1:length(m2dt.data)
    if m2dt.data(i,1)>0
        m2dt.data(i,1) = m2dt.data(i,1)-2*pi;
    end
end
m1dt.data = -m1dt.data;
m2dt.data = -m2dt.data;



%% Visualization
% MMCLab
Params.cboff = 1;
Params.Scale = -10;
Params.Th.P = 0; Params.Th.N = -Params.Th.P;
Params.Cmap = 'jet'; Params.PD = 1;
Params.EdgeColor = 'none';
PlotMeshSurface(m1mc,Params);  
view([90 0])
axis off
fig1 = gca;
saveas(fig1,[path,'/figures/Gmmc_intensity_',num2str(frequency),'MHz.png'])

Params1.cboff = 1;
Params1.Scale = 2*pi;
Params1.Th.P = 0; Params1.Th.N = -Params1.Th.P;
Params1.Cmap = 'hsv'; 
Params1.PD = 1;
Params1.OL = 0;
Params1.EdgeColor = 'none';
PlotMeshSurface(m2mc,Params1);   
view([90 0])
axis off
fig2 = gca;
saveas(fig2,[path,'/figures/Gmmc_phase_',num2str(frequency),'MHz.png'])


% NIRFASTer
Params.Scale = -10;
Params.Th.P = 0; Params.Th.N = -Params.Th.P;
Params.Cmap = 'jet'; Params.PD = 1;
Params.EdgeColor = 'none';
PlotMeshSurface(m1dt,Params);  
view([90 0])
axis off
fig3 = gca;
saveas(fig3,[path,'/figures/Gnirfaster_intensity_',num2str(frequency),'MHz.png'])

Params1.Scale = 2*pi;
Params1.Th.P = 0; Params1.Th.N = -Params1.Th.P;
Params1.Cmap = 'hsv'; 
Params1.PD = 1;
Params1.OL = 0;
Params1.EdgeColor = 'none';
PlotMeshSurface(m2dt,Params1);   
view([90 0])
axis off
fig4 = gca;
saveas(fig4,[path,'/figures/Gnirfaster_phase_',num2str(frequency),'MHz.png'])

