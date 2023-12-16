clear 
clc
close all

%% Load flux and initialize mcGs, mcGd
frequency = 1e8;     % Hz
path = '/data/eggebrecht/data1/Weihao/NeuroDOT_WF';
load([path,'/mesh.mat'])
load([path,'/MC/flux_wavelength1_source1.mat'])

data_s = zeros(size(mesh.nodes,1),size(flux.data,2),size(spos3,1),2);
data_d = zeros(size(mesh.nodes,1),size(flux.data,2),size(dpos3,1),2);
mcGs = zeros(2,size(spos3,1),size(mesh.nodes,1));
mcGd = zeros(2,size(dpos3,1),size(mesh.nodes,1));

for i = 1:size(spos3,1)
    load([path,'/MC/flux_wavelength1_source',num2str(i),'.mat'])
    data_s(:,:,i,1) = flux.data;
    load([path,'/MC/flux_wavelength2_source',num2str(i),'.mat'])
    data_s(:,:,i,2) = flux.data;
end

for i = 1:size(dpos3,1)
    load([path,'/MC/flux_wavelength1_detector',num2str(i),'.mat'])
    data_d(:,:,i,1) = flux.data;
    load([path,'/MC/flux_wavelength2_detector',num2str(i),'.mat'])
    data_d(:,:,i,2) = flux.data;
end

%% Make mcGs
for optode = 1:size(spos3,1)
    for wavelength=1:2
        data = data_s(:, :, optode, wavelength);
        if frequency == 0
            % convert time-resolved fluence to CW fluence
            fluence = sum(data*tstep,2); 
        else
            % convert time-resolved fluence to FD fluence
            w = 2*pi*frequency;
            [idx, it] = size(data);
            fluence = single(zeros(idx,1));
            rea = cos(w*(1:it)*tstep);
            ima = 1i*sin(w*(1:it)*tstep);
            for i=1:idx
                    fluence(i,1) = sum(data(i,1:it)*tstep.*rea(1:it)...
                              +data(i,1:it)*tstep.*ima(1:it));
            end
        end
        mcGs(wavelength,optode,:) = fluence(:,1); 
    end
end

%% Make mcGd
for optode = 1:size(dpos3,1)
    for wavelength=1:2
        data = data_d(:, :, optode, wavelength);
        if frequency == 0
            % convert time-resolved fluence to CW fluence
            fluence = sum(data*tstep,2); 
        else
            % convert time-resolved fluence to FD fluence
            w = 2*pi*frequency;
            [idx, it] = size(data);
            fluence = single(zeros(idx,1));
            rea = cos(w*(1:it)*tstep);
            ima = 1i*sin(w*(1:it)*tstep);
            for i=1:idx
                    fluence(i,1) = sum(data(i,1:it)*tstep.*rea(1:it)...
                              +data(i,1:it)*tstep.*ima(1:it));
            end
        end
        mcGd(wavelength,optode,:) = fluence(:,1); 
    end
end

save([path,'/MC/',num2str(frequency/1e6),'/GFunc_MC_',num2str(frequency/1e6),'MHz_.mat'], 'mcGs', 'mcGd', '-v7.3')

