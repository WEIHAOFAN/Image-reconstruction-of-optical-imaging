% clear
close all
clc

% frequency = 0;
% rs = 40;
% path = ['/data/eggebrecht/data1/Weihao/NeuroDOT_WF/'];
% load([path,'DT/',num2str(frequency),'/A_LUMO_36tiles_on_atlas_test.mat']);
% load([path,'mesh.mat'])
% load([path,'DT/',num2str(frequency),'/iA_lumo_0MHz.mat'])
% load('/data/eggebrecht/data1/matlab_codes/neurodot_dev/NeuroDOT_Beta/Support_Files/Spectroscopy/Ecoeff_Prahl.mat','prahlEC')
% keep1 = (info.pairs.WL==1 & info.pairs.r3d<=rs);
% keep2 = (info.pairs.WL==2 & info.pairs.r3d<=rs);
% Extco(1,:) = prahlEC(ismember(prahlEC(:,1),750),2:3).*[2.303,2.303]./10000000;
% Extco(2,:) = prahlEC(ismember(prahlEC(:,1),850),2:3).*[2.303,2.303]./10000000;

dim1 = 37;
dim2 = 13;
dim3 = 28;

surface_coor = mesh.nodes(mesh.bndvtx==1,:);
node_coor = change_space_coords([dim1, dim2, dim3],info.tissue.dim,'coord');
depth = min(pdist2(surface_coor, node_coor))

xx = info.tissue.dim.Good_Vox(30000,1);

pp = sub2ind([info.tissue.dim.nVx info.tissue.dim.nVy info.tissue.dim.nVz], dim1, dim2, dim3);
hemo_sim = zeros(2,size(a690,2));
hemo_sim(1,ismember(info.tissue.dim.Good_Vox, pp)) = 3.8;
hemo_sim(2,ismember(info.tissue.dim.Good_Vox, pp)) = -1.8;
rsd690 = squeeze(info.pairs.r3d(keep1));
rsd850 = squeeze(info.pairs.r3d(keep2));

rand1 = abs(randn(sum(keep1),1));
rand2 = abs(randn(sum(keep2),1)); 
[noise_real690,noise_img690,noise_real850,noise_img850] = realistic_noise(rsd690,rsd850,frequency,rand1,rand2);
if frequency==0
   noise690 = noise_real690;
   noise850 = noise_real850;
else
   noise690 = [noise_real690;noise_img690];
   noise850 = [noise_real850;noise_img850];
end
Xsim = 1e3*(Extco*hemo_sim)';
Ysim690_noise = a690*Xsim(:,1)+noise690;
Ysim850_noise = a850*Xsim(:,2)+noise850;
Xrecon690_noise = iA690*Ysim690_noise;
Xrecon850_noise = iA850*Ysim850_noise;
Xrecon_noise = [Xrecon690_noise';Xrecon850_noise'];
hemo_recon = Extco^(-1)*Xrecon_noise;

% Xsim = (Extco*hemo_sim)';
% Ysim690 = a690*Xsim(:,1);
% Ysim850 = a850*Xsim(:,2);
% Xrecon690 = iA690*Ysim690;
% Xrecon850 = iA850*Ysim850;
% Xrecon = [Xrecon690';Xrecon850'];
% hemo_recon = Extco^(-1)*Xrecon;

psf = Good_Vox2vol(hemo_recon(1,:),info.tissue.dim);
psf = psf./max(psf(:));
% pA_psf.slices_type='coord';
% pA_psf.PD = 1;pA_psf.Scale=1;pA_psf.Th.P=1e-2;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
% PlotSlices(t1,info.tissue.dim,pA_psf, psf)

psf = StrongestContiguousRegion1(psf,0.5);
psf(dim1,dim2,dim3) = 0.25;
[aFW, aFV, aLE, aER] = Cal_metrics_dev(dim1,dim2,dim3,psf);
recon = psf(psf>0);
    
pA_psf.slices_type='idx';pA_psf.slices=[dim1 dim2 dim3];pA_psf.view ='c';
pA_psf.PD = 1;pA_psf.Scale=1;pA_psf.Th.P=1e-2;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
PlotSlice(t1,info.tissue.dim,pA_psf, psf)
axis off
set(gcf,'position',[800 800 600 500])
fig1 = gca;
saveas(fig1,[path,'/figures/PSF_',num2str(floor(depth)),'_',num2str(frequency),'MHz_DTc.png'])

pA_psf.slices_type='idx';pA_psf.slices=[dim1 dim2 dim3];pA_psf.view ='t';
pA_psf.PD = 1;pA_psf.Scale=1;pA_psf.Th.P=1e-2;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
PlotSlice(t1,info.tissue.dim,pA_psf, psf)
axis off
set(gcf,'position',[800 800 600 700])
fig2 = gca;
saveas(fig2,[path,'/figures/PSF_',num2str(floor(depth)),'_',num2str(frequency),'MHz_DTt.png'])


% figure
% plot(1:length(Ysim690_noise), Ysim850_noise,'.-')
% 
% figure
% plot(rsd690, noise850,'.')

