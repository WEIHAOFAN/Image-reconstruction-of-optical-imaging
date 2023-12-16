%% Script for model generation given an initial segmented vokumetric mask and optode set

% Volumetric segmentation is currently not supported, but there are many
% options available (e.g., FreeSurfer). The starting point to this script
% is a segmented volume and a Pad file, both of which are contained within
% the Support_Files folder of the toolbox. To generate your own head model
% and light model, please just adapt the steps outlined below to your own
% volumetric mask and source-detector grid design.

clear 
close all
clc
% %% Load a Segmented Volume and generate head mesh
[mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111'],[],'4dfp');
p.Cmap='jet';p.Scale=5;p.Th.P=0;p.Th.N=-p.Th.P;p.PD=1;p.BG=[0,0,0];
PlotSlices(mask,infoT1,p)       % Visualize the segmented mask: 
%         % Note, PlotSlices is an interactive plot. 
%         % Hit q or the middle mouse button to quit
% 
% 
% % Parameters for generating your mesh
meshname='atlas';      % Provide a name for your mesh name here
% param.facet_distance=3;   % Node position error tolerance at boundary
% param.facet_size =1.2;      % boundary element size parameter
% param.cell_size=1.5;        % Volume element size parameter
% param.info=infoT1;
% param.Offset=[0,0,0];
% param.r0=5;                  % nodes outside of mask must be set to scalp==5;
% param.CheckMeshQuality=0;
% param.Mode=0;
% tic;mesh=NirfastMesh_Region(mask,meshname,param);toc
% 
pM.orientation='coord';pM.Cmap.P='gray';
% PlotMeshSurface(mesh,pM)
% %%% IF you get an error, this is due to NIRFAST using an old version of
% %%% mode. Go into NIRFAST/toolbox, change mode.m to modeOLD.m. Then the
% %%% correct version of mode will be used and this visualization will work.
% 
% % Put coordinates back in true space
% mesh.nodes=change_space_coords(mesh.nodes,infoT1,'coord');
% PlotMeshSurface(mesh,pM)
% 
% % Look at the inside
% m1.nodes=mesh.nodes;m1.elements=mesh.elements;
% m2=CutMesh(m1,find(m1.nodes(:,2)>0));
% [~,Ib]=ismember(m2.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
% m2.region=mesh.region(Ib);
% PlotMeshSurface(m2,pM);
% 
% 
% %% Load initial Pad file depending on data structure:
padname='LUMO_36tiles';
% load(['Pad_',padname,'.mat']); % Pad file
% tpos=cat(1,info.optodes.spos3,info.optodes.dpos3);
% Ns=size(info.optodes.spos3,1);
% Nd=size(info.optodes.dpos3,1);
% rad=info;
% PlotSD(tpos(1:Ns,:),tpos((Ns+1):end,:),'norm');
% PlotCap(info)
% 
% 
% %% View cap position with and witclear all
% 
% if ~exist('mesh','var')
% %     [mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111'],[],'4dfp');
%     mesh=load_mesh(meshname);
%     mesh.nodes=change_space_coords(mesh.nodes,infoT1,'coord');
%     PlotMeshSurface(mesh,pM)
% end
% 
% PlotMeshSurface(mesh,pM);PlotSD(tpos(1:Ns,:),tpos((Ns+1):end,:),'norm',gcf);
% 
% %% Move grid from arbitrary location to approximate target on mesh
% 
% % Adjust parameters by hand, update position, check
% tpos2=tpos;
% dx=0;
% dy=0;
% dz=-5;
% dS=1.2;
% dxTh=0;
% dyTh=0;
% dzTh=0;
% 
% tpos2=rotate_cap(tpos2,[dxTh,dyTh,dzTh]);
% tpos2(:,1)=tpos2(:,1)+dx;
% tpos2(:,2)=tpos2(:,2)+dy;
% tpos2(:,3)=tpos2(:,3)+dz;
% tpos2=scale_cap(tpos2,dS);
% 
pM.reg=0; pM.EdgeColor='none';
% paramsFoci.radius=2; paramsFoci.faces=10; paramsFoci.color=cat(1,repmat([1 0 0],Ns,1),repmat([0 0 1],Nd,1));
% PlotMeshSurface(mesh,pM);Draw_Foci_191203(tpos2,paramsFoci)
% 
% %% Relax grid on head
% m0.nodes=mesh.nodes;
% m0.elements=mesh.elements;
% spos3=tpos2(1:Ns,:);
% dpos3=tpos2((Ns+1):end,:);
% tposNew=gridspringfit_ND2(m0,rad,spos3,dpos3);
% 
% info.optodes.spos3=tposNew(1:Ns,:);
% info.optodes.dpos3=tposNew((Ns+1):end,:);
% info.tissue.infoT1=infoT1;
% info.tissue.affine=eye(4);
% info.tissue.affine_target='MNI';
% save(['Pad_',padname,'_on_',meshname,'.mat'],'info','tposNew'); % Pad file
% 
% 
% %% View optode positions
% PlotMeshSurface(mesh,pM);PlotSD(tposNew(1:Ns,:),tposNew((Ns+1):end,:),'render',gcf);
% PlotMeshSurface(mesh,pM);PlotSD(spos3,dpos3,'render',gcf);
% 
% %% PREPARE! --> mesh and grid in coord space
gridname=padname;
% mesh=PrepareMeshForNIRFAST(mesh,[meshname,'_',gridname],tposNew);
% PlotMeshSurface(mesh,pM);PlotSD(mesh.source.coord(1:Ns,:),...
%     mesh.source.coord((Ns+1):end,:),'render',gcf);
% view([90 0])
% 
% % One last visualization check...
% m3=CutMesh(mesh,intersect(find(mesh.nodes(:,3)>0),find(mesh.nodes(:,1)>0)));
% [Ia,Ib]=ismember(m3.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
% m3.region=mesh.region(Ib);
% PlotMeshSurface(m3,pM);PlotSD(mesh.source.coord(1:Ns,:),...
%     mesh.source.coord((1+Ns):end,:),'render',gcf);
% 
% spos3=tposNew(1:Ns,:);
% dpos3=tposNew((Ns+1):end,:);
% save('mesh.mat','mesh','spos3','dpos3','-v7.3')
load('/data/eggebrecht/data1/Weihao/NeuroDOT_WF/mesh.mat')
load('/data/eggebrecht/data1/Weihao/NeuroDOT_WF/Pad_LUMO_36tiles_on_atlas.mat')
Ns = size(spos3,1);
Nd = size(dpos3,1);
PlotMeshSurface(mesh,pM);PlotSD(spos3,dpos3,'render',gcf);

%% Calculate Sensitivity Profile
% Parameters
flags.tag=[padname,'_on_',meshname,'_test'];
flags.gridname=gridname;
flags.meshname=meshname;
flags.head='info';
flags.info=infoT1;                  % Your T1 info file
flags.gthresh=1e-3;                 % Voxelation threshold in G
flags.voxmm=2;                      % Voxelation resolution (mm)
flags.labels.r1='csf';             % Regions for optical properties
flags.labels.r2='white';
flags.labels.r3='gray';
flags.labels.r4='bone';
flags.labels.r5='skin';
flags.op.lambda=[750,850];          % Wavelengths (nm)
flags.op.mua_skin=[0.0183,0.0190];  % Baseline absorption
flags.op.mua_bone=[0.0118,0.0139];
flags.op.mua_csf=[0.0040,0.0040];
flags.op.mua_gray=[0.0201,0.0192];
flags.op.mua_white=[0.0171,0.0208];
flags.op.musp_skin=[0.8,0.64];     % Baseline reduced scattering coeff
flags.op.musp_bone=[1.0,0.84];
flags.op.musp_csf=[0.3,0.3];
flags.op.musp_gray=[0.9727,0.6726];
flags.op.musp_white=[1.3333,1.0107];
flags.op.n_skin=[1.4,1.4];          % Index of refration
flags.op.n_bone=[1.4,1.4];
flags.op.n_csf=[1.4,1.4];
flags.op.n_gray=[1.4,1.4];
flags.op.n_white=[1.4,1.4];
flags.srcnum=Ns;                    % Number of sources
flags.t4=eye(4);                    % T1/dim to MNI atlas *** change this to register your vol to atlas
flags.t4_target='MNI'; % string
flags.makeA=0; % don't make A, just make G
flags.Hz=4e8;
if flags.Hz, flags.tag = [flags.tag,'FD']; end

Ti=tic;[~,dim,~,Gs,Gd,dc]=makeAnirfaster(mesh,flags); % size(A)= [Nwl, Nmeas, Nvox]
disp(['<makeAnirfast took ',num2str(toc(Ti))])


%% Package data and save A with only 1st through 5th nearest neighbors to save space

% make A with measurements <5cm
keep=info.pairs.r3d<=50;
temp=struct;
temp.Src=info.pairs.Src(keep);
temp.Det=info.pairs.Det(keep);
temp.NN=info.pairs.NN(keep);
temp.WL=info.pairs.WL(keep);
temp.lambda=info.pairs.lambda(keep);
temp.Mod=info.pairs.Mod(keep);
temp.r2d=info.pairs.r2d(keep);
temp.r3d=info.pairs.r3d(keep);
info.pairs=temp;
flags.infoA=info;

[A,Gsd]=g2a_200225(Gs,Gd,dc,dim,flags);

info.tissue.dim=dim;
info.tissue.affine=flags.t4;
info.tissue.infoT1=infoT1;
info.tissue.affine_target='MNI';
info.tissue.flags=flags;

info.tissue.dim.center = info.tissue.dim.center - [-2, 0, 0];
t1=affine3d_img(mask,infoT1,info.tissue.dim,eye(4)); % put anatomical volume in dim space
% t = toc(t);

save(['A_',flags.tag,'.mat'],'A','info','t1','Gsd','-v7.3')



% % Visualize aspects of sensitivity profile
% 
% 
% keep=info.pairs.WL==2 & info.pairs.Src==78 & info.pairs.Det==80;
% foo=squeeze(A(keep,:));              % Single meas pair
% fooV=Good_Vox2vol(foo,info.tissue.dim);
% fooV=fooV./max(fooV(:));
% fooV=log10(1e3.*fooV);                  % top 2 o.o.m.
% pA.PD=1;pA.Scale=3;pA.Th.P=0;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,info.tissue.dim,pA,fooV)
% 
% % FFR
% keep=(info.pairs.WL==2 & info.pairs.NN<=4);
% a=squeeze(A(keep,:));
% iA=Tikhonov_invert_Amat(a,0.01,0.1);
% iA=smooth_Amat(iA,info.tissue.dim,3);
% ffr=makeFlatFieldRecon(a,iA);
% % 
% fooV=Good_Vox2vol(ffr,info.tissue.dim);
% fooV=fooV./max(fooV(:));
% pA.PD=1;pA.Scale=1;pA.Th.P=1e-4;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,info.tissue.dim,pA,fooV)

% %SurfPlot
% [t4,scale,header] = loadAviT4([],'TC35617_dim_to_MNI_t4');
% ffrMNI = affine3d_img(fooV,dim,infoT1,t4);
% pA.view = 'post'; pA.ctx = 'std'; % pA.view = 'post', 'lat', 'dorsal'; pA.ctx = 'inf'; 'vinf'
% load('MNI164k_big.mat')
% PlotInterpSurfMesh(ffrMNI,MNIl,MNIr,infoT1,pA);
%  
% load('A_AdultV24x28_on_Example_Mesh_weihao.mat','A','info','t1')
% frequency = 0;
% keep = (info.pairs.WL==2 & info.pairs.NN<5);
% a = squeeze(A(keep,:));

