function [output_norm, output_signal, output_noise, output_snr] = Cal_unthreshed_metrics_neuroDOT(x_perturb, y_perturb, z_perturb, PSF_raw,ff, pp,  N_recon)
%% calculate image norm, signal, noise and SNR
[dim1, dim2, dim3] = size(PSF_raw);
ax = ceil(dim1/2);
ay = ceil(dim2/2);
perturb = zeros(dim1, dim2, dim3);
perturb(x_perturb, y_perturb, z_perturb) = pp;
psf = StrongestContiguousRegion1(PSF_raw,0.5);
matrixn = pp*PSF_raw./max(PSF_raw(:));
error = reshape((matrixn-perturb),dim1*dim2*dim3,1);
anorm = norm(error);
asignal = mean(PSF_raw(psf>0));
imagenoise = N_recon.*ff;
anoise = std(imagenoise(:));
asnr = asignal/anoise;

output_norm = anorm;
output_signal = asignal;
output_noise = anoise;
output_snr = asnr;

end