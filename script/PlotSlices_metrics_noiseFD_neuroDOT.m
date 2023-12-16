function [Output_FW, Output_FV, Output_LE, Output_ER, Output_Norm, Output_Signal, Output_Noise, Output_SNR] = PlotSlices_metrics_noiseFD_neuroDOT(fooV, info, keep1, keep2, aWL1, aWL2, iAWL1, iAWL2, frequency, ROI, Extco, p_HbO, p_HbR)
%% Plot the FWHM and Localization error in 3 slices. Using the sensitivity
%%matrices to get the Xrecon in 3D space, which is fooV_psf. Then get the
%%metrics, FWHM and Localization error. Finally plot them onto underlay mesh. 

%% Calculate the metrics i different voxel.
p_HbT = p_HbO+p_HbR;
[fooV_x, fooV_y, fooV_z] = size(fooV);
l_FW = zeros(fooV_x, fooV_y, fooV_z, 3);
l_FV = zeros(fooV_x, fooV_y, fooV_z, 3);
l_LE = zeros(fooV_x, fooV_y, fooV_z, 3);
l_ER = zeros(fooV_x, fooV_y, fooV_z, 3);
l_Norm = zeros(fooV_x, fooV_y, fooV_z, 3);
l_Signal = zeros(fooV_x, fooV_y, fooV_z, 3);
l_Noise = zeros(fooV_x, fooV_y, fooV_z, 3);
l_SNR = zeros(fooV_x, fooV_y, fooV_z, 3);
for m=1:size(info.tissue.dim.Good_Vox,1)
    [i,j,k] = ind2sub([fooV_x, fooV_y, fooV_z],info.tissue.dim.Good_Vox(m,1));
    [hemo_recon_noise, image_noise] = MakeFD_HemoreconfooV_noise_neuroDOT(info.tissue.dim.Good_Vox(m,1),info,keep1,keep2,aWL1,aWL2,iAWL1,iAWL2,frequency,Extco, p_HbO, p_HbR);     
    HbO_psf=Good_Vox2vol(hemo_recon_noise(1,:),info.tissue.dim);
    HbR_psf=Good_Vox2vol(hemo_recon_noise(2,:),info.tissue.dim);
    HbT_psf=HbO_psf+HbR_psf;
    HbO_noise=Good_Vox2vol(image_noise(1,:),info.tissue.dim);
    HbR_noise=Good_Vox2vol(image_noise(2,:),info.tissue.dim);
    HbT_noise=HbO_noise+HbR_noise;
    [l_Norm(i,j,k,1), l_Signal(i,j,k,1),l_Noise(i,j,k,1), l_SNR(i,j,k,1)] = Cal_unthreshed_metrics_neuroDOT(i,j,k,HbO_psf,ROI,p_HbO,HbO_noise);
    [l_Norm(i,j,k,2), l_Signal(i,j,k,2),l_Noise(i,j,k,2), l_SNR(i,j,k,2)] = Cal_unthreshed_metrics_neuroDOT(i,j,k,HbR_psf,ROI,p_HbR,HbR_noise);
    [l_Norm(i,j,k,3), l_Signal(i,j,k,3),l_Noise(i,j,k,3), l_SNR(i,j,k,3)] = Cal_unthreshed_metrics_neuroDOT(i,j,k,HbT_psf,ROI,p_HbT,HbT_noise);               
    HbO_psf=HbO_psf./abs(max(abs(HbO_psf(:))));
    HbR_psf=-HbR_psf./abs(max(abs(HbR_psf(:))));
    HbT_psf=HbT_psf./abs(max(abs(HbT_psf(:))));
    HbO_psf = StrongestContiguousRegion(HbO_psf,0.5);
    HbR_psf = StrongestContiguousRegion(HbR_psf,0.5);
    HbT_psf = StrongestContiguousRegion(HbT_psf,0.5);
    [l_FW(i,j,k,1), l_FV(i,j,k,1), l_LE(i,j,k,1), l_ER(i,j,k,1)] = Cal_metrics_dev(i,j,k,HbO_psf);
    [l_FW(i,j,k,2), l_FV(i,j,k,2), l_LE(i,j,k,2), l_ER(i,j,k,2)] = Cal_metrics_dev(i,j,k,HbR_psf);
    [l_FW(i,j,k,3), l_FV(i,j,k,3), l_LE(i,j,k,3), l_ER(i,j,k,3)] = Cal_metrics_dev(i,j,k,HbT_psf);
    if mod(m,1e3)==0
        m
    end
end

Output_FW = l_FW;
Output_FV = l_FV;
Output_LE = l_LE;
Output_ER = l_ER;
Output_Norm = l_Norm;
Output_Signal = l_Signal;
Output_Noise = l_Noise;
Output_SNR = l_SNR;
end