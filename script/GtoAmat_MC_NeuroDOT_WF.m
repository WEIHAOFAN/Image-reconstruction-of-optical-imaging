function [A,dim,Gsd,mcGs,mcGd,dc]=GtoAmat_MC_NeuroDOT_WF(mcGs,mcGd,mesh,dc,flags,dtGs,dtGd)

% This function take a set of Green's functions from MMCLab and NIRFAST and a mesh and creates an
% A-matrix.
%
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.

%% Paramters
flags.keepmeth='glevel';
if ~isfield(flags,'gthresh'),flags.gthresh=10^-3;end
if ~isfield(flags,'voxmm'),flags.voxmm=2;end
if ~isfield(flags,'makeA'),flags.makeA=1;end
flags.GV=1;


%% Interpolate  
disp('>Finding Voxel Grid and Cropping Limits')
[vox,dim]=getvox(mesh.nodes,cat(2,dtGs,dtGd),flags);

disp('>Finding Mesh to Voxel Converstion')
tic;TR=triangulation(mesh.elements,mesh.nodes);toc
tic;[t,p] = pointLocation(TR,reshape(vox,[],3));toc

% mytsearchn takes too long
% tic;[t,p] = mytsearchn(mesh,reshape(vox,[],3));toc
% tsearchn is unstable if possible nonconvexity...  : /
% [t,p] = tsearchn(mesh.nodes,mesh.elements,reshape(vox,[],3)); 

disp('>Interpolating Greens Functions')
mcGs=voxel(mcGs,t,p,mesh.elements);
mcGd=voxel(mcGd,t,p,mesh.elements);

disp('>Interpolating Optical Properties')
dc=voxel(dc,t,p,mesh.elements);


%% Create dim.Good_Vox
if flags.GV==1
    [~,~,dc,dim]=Make_Good_Vox(dtGs,dtGd,dc,dim,mesh,flags);
end

% %% Save voxellated G etc
% disp('Saving Voxellated Green''s Functions')
save(['GFunc_',flags.tag,'_VOX.mat'],'mcGs','mcGd','dim','flags',...
    't','p','dc','vox','-v7.3')
% clear vox p t mesh


%% Create A-matrix
if flags.makeA
    disp('>Making A-Matrix')
    [A,Gsd]=g2a(mcGs,mcGd,dc,dim,flags);
    disp('>Saving A-Matrix')
    save(['A_',flags.tag,'.mat'],'A','dim','flags','Gsd','-v7.3')
    
else
    A=[];
    Gsd=[];
end