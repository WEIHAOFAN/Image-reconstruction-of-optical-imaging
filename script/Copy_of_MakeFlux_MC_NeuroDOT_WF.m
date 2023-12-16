clear 
close all
clc

% The script simulate the flux of each optode with two wavelengths using MMCLAB.
% The script require the mesh and pad files right after PrepareforNirfast.

%% load mask, mesh, and pad files
[mask,infoT1]=LoadVolumetricData(['Segmented_MNI152nl_on_MNI111'],[],'4dfp');        
load('/data/eggebrecht/data1/Weihao/NeuroDOT_WF/mesh.mat')
load('/data/eggebrecht/data1/Weihao/NeuroDOT_WF/Pad_LUMO_36tiles_on_atlas.mat')

%% Initialization
meshname='atlas';         %  name of the mesh
padname='LUMO_36tiles';        %  name of the pad
gridname=padname;

% Parameters prepare for MMCLAB
elemprop = element_region(mesh.region,mesh.elements);
cfg.nphoton=5e8;    % photon numbers launched for each optode
cfg.node=double(mesh.nodes);
cfg.elem=double(mesh.elements);
cfg.elemprop=double(elemprop);

% time windows
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=1e-10;
cfg.gpuid=1;

% other simulation parameters
cfg.isreflect=1;
cfg.debuglevel='TP';
cfg.method='elem';      %  MMC mode
cfg.srctype='isotropic';


%% Check the mesh and pad
pM.orientation='coord'; pM.Cmap.P='gray'; pM.reg=0; pM.EdgeColor='none';
PlotMeshSurface(mesh,pM);PlotSD(spos3,dpos3,'render',gcf);

tic
%%  Calculate Sensitivity Profile by MMC
for wavelength = 1
    for i = 232:size(tposNew,1)
        clear data;
        clear flux;

        % medium labels: 0-ambient air, 1-air cavities, 2-scalp, 3-skull,
        %                          4-csf, 5-gray matter, 6-white matter
        if wavelength == 1
            cfg.prop=[0,0,1,1;               %  air
                    0.0040 0.01 0.9 1.4;            %  CSF
                    0.0171 13.333 0.9 1.4;       %  white matter
                    0.0201 9.727 0.9 1.4;        %  gray matter
                    0.0118 10.0 0.9 1.4          %  bone
                    0.0183 8.0 0.9 1.4];          %  skin
                
        else
            cfg.prop=[0,0,1,1;               %  air
                    0.0040 0.01 0.9 1.4;            %  CSF
                    0.0208 10.107 0.9 1.4;       %  white matter
                    0.0192 6.726 0.9 1.4;        %  gray matter
                    0.0139 8.4 0.9 1.4          %  bone
                    0.0190 6.4 0.9 1.4];          %  skin
        end
  
        % light direction
        surface = mesh.nodes(mesh.bndvtx==1,:);
        cfg.srcpos = mesh.source.coord(i,:);
        [nx,ny,nz] = headnorm(surface, mesh.source.coord(i,:));    % calculate the norm vector
        cfg.srcdir = [nx,ny,nz]; %inward-pointing source
        mesh.source.coord(i,:) = mesh.source.coord(i,:)+0.01*cfg.srcdir;
        
        % model the Green's functions (flux.data)
        [flux]=mmclab(cfg);
        tstep = cfg.tstep;
        
        % save the flux
        if i<=size(info.optodes.spos2,1)
            save(['flux_wavelength',num2str(wavelength),'_source',num2str(i),'.mat'],'flux','tstep','-v7.3')
        else
            save(['flux_wavelength',num2str(wavelength),'_detector',num2str(i-size(info.optodes.spos2,1)),'.mat'],'flux','tstep','-v7.3')
        end

    end
end
toc
