%% The script simulate the image reconstruction based on the modelde activiation 

% The script uses given sensitivity matrix, mesh and extinction coefficient matrix  
% to generate simulated light data, y, and further reconstruct the brain functions. 
% The simulated reconstruction is then processed for calculation of image quality 
% metrics changing with depth.

clear 
close all
clc

%% Initialization
frequency = 300;
hemo = 1;   % selct the hemoglobin types  (1)HbO (2)HbR (3)HbT
noise = 1;   % 0 for noise free
l1 = 0.01;   % Tikhonov regularization
l2 = 0.1;    % Spatially variant regularization
smooth_p = 5;     % smooth parameter 
p_HbO = 3.8;    % point perturbation of HbO
p_HbR = -1.8;   % point perturbation of HbR
thresh = 7;    % threshold parameter of ROI

%% load A matrix, mesh, and extinction coefficient matrix
path = '/data/eggebrecht/data1/Weihao/NeuroDOT_WF/';
load([path,'DT/300/A_LUMO_36tiles_on_atlas_testFD.mat'])
% load([path,'MC/0/A_MC_LUMO_36tiles_on_atlas_test.mat'])
load([path,'/mesh.mat'])
load([path,'/Pad_LUMO_36tiles_on_atlas.mat'],'tposNew')
load('/data/eggebrecht/data1/matlab_codes/neurodot_dev/NeuroDOT_Beta/Support_Files/Spectroscopy/Ecoeff_Prahl.mat','prahlEC')

% set extinction coefficient
WL = unique(info.pairs.lambda);
WLmod = mod(WL,2);
Extco(find(WLmod==1),:) = (prahlEC(ismember(prahlEC(:,1),(WL(WLmod==1)-1)),2:3)+...
           prahlEC(ismember(prahlEC(:,1),(WL(WLmod==1)+1)),2:3)).*[2.303,2.303]./2e7;
Extco(find(WLmod==0),:) = prahlEC(ismember(prahlEC(:,1),WL(WLmod==0)),2:3).*[2.303,2.303]./1e7;
% the 1e7/2e7 factors are intended to convert the unit to be 1/(mm*uM) 

% select measurements and wavelength
keep1 = (info.pairs.WL==1 & info.pairs.r3d<=40);
keep2 = (info.pairs.WL==2 & info.pairs.r3d<=40);

%% Invert A
aWL1 = squeeze(A(keep1,:));
aWL2 = squeeze(A(keep2,:));
if frequency   % FD mode
    aWL1 = [real(aWL1); imag(aWL1)];
    aWL2 = [real(aWL2); imag(aWL2)];
end

% invert A matrix
iAWL1=Tikhonov_invert_Amat(aWL1,l1,l2);
iAWL2=Tikhonov_invert_Amat(aWL2,l1,l2);
iAWL1=smooth_Amat(iAWL1,info.tissue.dim,smooth_p);
iAWL2=smooth_Amat(iAWL2,info.tissue.dim,smooth_p);

%% check flat field region (FFR)
ffrWL1=makeFlatFieldRecon(aWL1,iAWL1);
ffrWL2=makeFlatFieldRecon(aWL2,iAWL2);

% WL1 for example
fooV=Good_Vox2vol(ffrWL1,info.tissue.dim);
fooV=fooV./max(fooV(:));
pA.PD=1;pA.Scale=1;pA.Th.P=1e-2;pA.Th.N=-pA.Th.P;pA.CH=0;
% PlotSlices(t1,info.tissue.dim,pA,fooV)

% WL2 for example
fooV=Good_Vox2vol(ffrWL2,info.tissue.dim);
fooV=fooV./max(fooV(:));
pA.PD=1;pA.Scale=1;pA.Th.P=1e-2;pA.Th.N=-pA.Th.P;pA.CH=0;
% PlotSlices(t1,info.tissue.dim,pA,fooV)

%% save iA  (iA could be large. Don't save unless neccessary)
% save(['iA_',num2str(frequency),'Hz.mat'],'iAWL1','iAWL1','t1','aWL1','aWL2','fooV','-v7.3')

%% visualize medians of each metric changing with depth
% kick the first outer layers of optodes
cap_coor = mesh.source.coord;
outsource = [78, 79, 84, 85, 90, 89, 62, 55, 76, 1, 2, 5, 8, 11, 13, 14, 36, 28,...
             48, 49, 54, 37, 39, 40, 45, 46, 92, 95, 98, 100, 101, 104, 107, 91]; 
outdetector = size(info.optodes.spos2,1) + [117, 86, 82, 77, 74, 97, 104, 2, 18,... 
                          22, 48, 44, 40, 33,64, 69, 52, 57, 122, 129, 134, 141];
outcap = [outsource, outdetector]';
cap_coor(outcap(:),:) = [];

% caculate the node depths
surface_coor = mesh.nodes(mesh.bndvtx==1,:);
[x, y, z] = ind2sub([info.tissue.dim.nVx,info.tissue.dim.nVy,info.tissue.dim.nVz],info.tissue.dim.Good_Vox');
node_coor = change_space_coords([x; y; z]',info.tissue.dim,'coord');
depth = min(pdist2(surface_coor, node_coor));
% caculate the minimum distance of each node to every optodes
optd = min(pdist2(cap_coor, node_coor));

% Defeine ROI by cutting off boundary with optd-depth > threshold
roi = ones(size(node_coor,1),1);
roi((optd-depth)>thresh) = 0;
ROI = Good_Vox2vol(roi,info.tissue.dim);

%% Simulate PSFs within Good_vox and calculate image quality metrics
if noise && frequency
    [FW, FV, LE, ER, Norm, Signal, Noise, SNR] = PlotSlices_metrics_noiseFD_neuroDOT(fooV, info, keep1, keep2, aWL1, aWL2, iAWL1, iAWL2, frequency, ROI, Extco, p_HbO, p_HbR);
end
if noise && ~frequency
    [FW, FV, LE, ER, Norm, Signal, Noise, SNR] = PlotSlices_metrics_noise_neuroDOT(fooV, info, keep1, keep2, aWL1, aWL2, iAWL1, iAWL2, frequency, ROI, Extco, p_HbO, p_HbR);
end
if ~noise 
    [FW, FV, LE, ER, Norm, Signal, Noise, SNR] = PlotSlices_metrics_neuroDOT(fooV, info, aWL1, aWL2, iAWL1, iAWL2, frequency, ROI, Extco, p_HbO, p_HbR);
end
    
% save([path,'/MC/',num2str(frequency),'/metrics_with_',num2str(frequency),'Hz.mat'],'FV','FW','LE','ER','ROI','roi','depth')
save([path,'/DT/',num2str(frequency),'/metrics_with_',num2str(frequency),'Hz.mat'],'FV','FW','LE','ER','ROI','roi','depth')
 
% Volumetric views of image quality metrics
pA_psf.slices_type='coord';
pA_psf.PD = 1;pA_psf.Scale=1;pA_psf.Th.P=1e-2;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
PlotSlices(t1,info.tissue.dim,pA_psf, LE(:,:,:,hemo).*ROI)

%% Visualize image quality metrics
LE = LE(:,:,:,hemo)*info.tissue.flags.voxmm.*ROI;
FW = FW(:,:,:,hemo)*info.tissue.flags.voxmm.*ROI;
FV = FV(:,:,:,hemo)*info.tissue.flags.voxmm.*ROI;
ER = ER(:,:,:,hemo)*info.tissue.flags.voxmm.*ROI;
LE0 = LE(LE>0);
FW0 = FW(LE>0 & LE<=8);
FV0 = FV(LE>0 & LE<=8);
ER0 = ER(LE>0 & LE<=8);
depth1 = depth(roi>0);      % depth for LE and SR
depth2 = depth1(LE0<=8);    % depth for FWHM, FVHM and ER

% localization error
medianLE = zeros(35,1);
for i=1:1:35
    medianLE(i,1) = median(LE0(depth1(:)<i & depth1(:)>=(i-1)));
end
fig1 = figure;
plot(0.5:1:34.5,medianLE,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
xlim([0 25])
ylim([0 25])
ylabel('localization error (mm)')
xlabel('depth (mm)')
grid on
set(gca,'FontName','Arial','fontsize',14,'LineWidth',2,'gridlinestyle','--')

% success rate
SR = zeros(70,1);
for i=1:1:70
    SR(i,1) = 100*sum(LE0(depth1(:)<i/2 & depth1(:)>=(i-1)/2)<=8)/sum(LE0(depth1(:)<i/2 & depth1(:)>=(i-1)/2)>0);
end
fig2 = figure;
plot(0.25:0.5:34.75,smooth(SR),'.-','MarkerSize',8,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
xlim([0 30])
ylim([0 100])
ylabel('success rate %')
xlabel('depth (mm)')
grid on
set(gca,'FontName','Arial','fontsize',14,'LineWidth',2,'gridlinestyle','--')

% FWHM
medianFW = zeros(35,1);
for i=1:1:35
    medianFW(i,1) = median(FW0(depth2(:)<i & depth2(:)>=(i-1)));
end
fig3 = figure;
plot(0.5:1:34.5,medianFW,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
xlim([0 25])
% ylim([0 20])
ylabel('FWHM (mm)')
xlabel('depth (mm)')
grid on
set(gca,'FontName','Arial','fontsize',14,'LineWidth',2,'gridlinestyle','--')

% FVHM
medianFV = zeros(35,1);
for i=1:1:35
    medianFV(i,1) = median(FV0(depth2(:)<i & depth2(:)>=(i-1)));
end
fig4 = figure;
plot(0.5:1:34.5,medianFV,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
% xlim([0 25])
% ylim([0 20])
ylabel('FVHM (mm)')
xlabel('depth (mm)')
grid on
set(gca,'FontName','Arial','fontsize',14,'LineWidth',2,'gridlinestyle','--')

% ER
medianER = zeros(35,1);
for i=1:1:35
    medianER(i,1) = median(ER0(depth2(:)<i & depth2(:)>=(i-1)));
end
fig5 = figure;
plot(0.5:1:34.5,medianER,'o-','MarkerSize',6,'MarkerEdgeColor',[0.9,0.2,0.2],...
            'MarkerFaceColor',[0.9,0.2,0.2],'LineWidth',2,'Color',[0.9,0.2,0.2])
xlim([0 25])
% ylim([0 20])
ylabel('Effective resolution (mm)')
xlabel('depth (mm)')
grid on
set(gca,'FontName','Arial','fontsize',14,'LineWidth',2,'gridlinestyle','--')

% save(['Depth chanegs of metrics_with_',num2str(frequency),'Hz.mat'],'depth1','depth2','medianLE','SR','medianFW','medianFV',,'medianER','-v7.3')

