clear
clc
close all


%% Load mesh, Green's function by NIRFAST and MMCLab, and mask
frequency = 0;
path = '/data/eggebrecht/data1/Weihao/NeuroDOT_WF';
load([path,'/mesh.mat'],'mesh');
load([path,'/Pad_LUMO_36tiles_on_atlas.mat'],'info');
load([path,'/MC/',num2str(frequency),'/GFunc_MC_',num2str(frequency/1e6),'MHz_.mat'])   % MMCLab
[mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111'],[],'4dfp');

%% Initialization
srcnum = size(mcGs,1);
detnum = size(mcGd,1);
numpt = size(mcGs,3);
measnum = srcnum*detnum;
meshname='atlas'; 
padname='LUMO_36tiles';
gridname=padname;

%% Settings for GtoA
% Flags parameters
flags.tag=[padname,'_on_',meshname,'_test'];
flags.gridname=gridname;
flags.meshname=meshname;
flags.head='info';
flags.info=infoT1;                  % Your T1 info file
flags.gthresh=1e-3;                 % Voxelation threshold in G
flags.voxmm=2;                      % Voxelation resolution (mm)
flags.op.lambda=[750,850];          % Wavelengths (nm)
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
flags.srcnum=srcnum;                    % Number of sources
flags.t4=eye(4);                    % T1/dim to MNI atlas *** change this to register your vol to atlas
flags.t4_target='MNI'; % string
flags.makeA=0; % don't make A, just make G
flags.Hz=frequency;


%% Parameters initialization
labels=flags.labels;
op=flags.op;
numc=length(op.lambda);

f=fieldnames(labels);
region=mesh.region;
regtypes=unique(region);

% initialize arrays
a=zeros(size(region));
kappa=zeros(size(region));
ri=zeros(size(region));
numNodes=size(mesh.nodes,1);
dc=zeros(numc,numNodes);


%%  Choose wavelengths and get Gfuncts then save them
for lambda=1:2
    for n=regtypes'
        roi=find(region==n);
        
        floop=['r',num2str(n)];
        
        if ~isfield(labels,floop)
            error(['** There is not description matching region number ',...
                num2str(n),' **'])
        elseif ~isfield(op,['mua_',labels.(floop)]) || ...
                ~isfield(op,['musp_',labels.(floop)]) || ...
                ~isfield(op,['n_',labels.(floop)])
            error(['** Not all optical properties are specified for region ',...
                labels.(floop),' **'])
        end
        
        a(roi)=op.(['mua_',labels.(floop)])(lambda);                % mua
        kappa(roi)=1/(3*(op.(['mua_',labels.(floop)])(lambda)+...   % D.C.
            op.(['musp_',labels.(floop)])(lambda)));
        ri(roi)=op.(['n_',labels.(floop)])(lambda);                 % Ind. of Ref.
    end
    
    mesh.mua=a;
    mesh.kappa=kappa;
    mesh.ri=ri;

    dc(lambda,:)=mesh.kappa;
end

[~,dim,~,mcGs,mcGd,dc]= GtoAmat(mcGs,mcGd,mesh,dc,flags);
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
[A,mcGsd]=g2a_200225(mcGs,mcGd,dc,dim,flags);

info.tissue.dim=dim;
info.tissue.affine=flags.t4;
info.tissue.infoT1=infoT1;
info.tissue.affine_target='MNI';
info.tissue.flags=flags;
info.tissue.dim.center = info.tissue.dim.center - [-2, 0, 0];
t1=affine3d_img(mask,infoT1,info.tissue.dim,eye(4)); % put anatomical volume in dim space

% keep=info.pairs.WL==2 & info.pairs.Src==1 & info.pairs.Det==1;
% foo=squeeze(A(keep,:));              % Single meas pair
% fooV=Good_Vox2vol(foo,dim);
% fooV=fooV./max(fooV(:));
% fooV=log10(1e3.*fooV);                  % top 2 o.o.m.
% pA.PD=1;pA.Scale=3;pA.Th.P=0;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,dim,pA,fooV)

% Save A matrix
save([path,'/MC/',num2str(frequency),'/A_MC_',flags.tag,'.mat'],'A','info','t1','mcGsd','-v7.3')


