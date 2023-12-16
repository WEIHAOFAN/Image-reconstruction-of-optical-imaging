function [Output_FW, Output_FV, Output_LE, Output_ES] = PlotSlices_metrics_neuroDOT(fooV, info, a690, a850, iA690, iA850, Extco, p_HbO, p_HbR)
%% Plot the FWHM and Localization error in 3 slices. Using the sensitivity
%%matrices to get the Xrecon in 3D space, which is fooV_psf. Then get the
%%metrics, FWHM and Localization error. Finally plot them onto underlay mesh. 

%% Calculate the metrics i different voxel.
[fooV_x, fooV_y, fooV_z] = size(fooV);
l_FW = zeros(fooV_x, fooV_y, fooV_z, 3);
l_FV = zeros(fooV_x, fooV_y, fooV_z, 3);
l_LE = zeros(fooV_x, fooV_y, fooV_z, 3);
l_ES = zeros(fooV_x, fooV_y, fooV_z, 3);
% for m=1:size(info.tissue.dim.Good_Vox,1)
for m=20000
    [i,j,k] = ind2sub([fooV_x, fooV_y, fooV_z],info.tissue.dim.Good_Vox(m,1));
    hemo_recon = MakeFD_HemoreconfooV_neuroDOT(info.tissue.dim.Good_Vox(m,1),info,a690,a850,iA690,iA850,Extco, p_HbO, p_HbR);
    HbO_psf=Good_Vox2vol(hemo_recon(1,:),info.tissue.dim);
    HbR_psf=Good_Vox2vol(hemo_recon(2,:),info.tissue.dim);
    HbT_psf=HbO_psf+HbR_psf;
    HbO_psf=HbO_psf./abs(max(abs(HbO_psf(:))));
    HbR_psf=-HbR_psf./abs(max(abs(HbR_psf(:))));
    HbT_psf=HbT_psf./abs(max(abs(HbT_psf(:))));
    HbO_psf = StrongestContiguousRegion1(HbO_psf,0.5);
    HbR_psf = StrongestContiguousRegion1(HbR_psf,0.5);
    HbT_psf = StrongestContiguousRegion1(HbT_psf,0.5);
    [l_FW(i,j,k,1), l_FV(i,j,k,1), l_LE(i,j,k,1), l_ES(i,j,k,1)] = Cal_metrics_dev(i,j,k,HbO_psf);
    [l_FW(i,j,k,2), l_FV(i,j,k,2), l_LE(i,j,k,2), l_ES(i,j,k,1)] = Cal_metrics_dev(i,j,k,HbR_psf);
    [l_FW(i,j,k,3), l_FV(i,j,k,3), l_LE(i,j,k,3), l_ES(i,j,k,1)] = Cal_metrics_dev(i,j,k,HbT_psf);
end

Output_FW = l_FW;
Output_FV = l_FV;
Output_LE = l_LE;
Output_ES = l_ES;

end