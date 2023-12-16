function [hemo_recon_noise, image_noise] = MakeFD_HemoreconfooV_noise_neuroDOT(perturb, info, keep1, keep2, a690, a850, iA690, iA850, frequency, Extco, p_HbO, p_HbR)
%% Input the perturbation coordinates (dim1 dim2 dim3), and sensitivity
%%matrcies A and iA_psf, return the output matrix, which is the Xrecon in
%%3D space.

%Make noise perturbation
hemo_sim = zeros(2,size(a690,2));
hemo_sim(1,ismember(info.tissue.dim.Good_Vox, perturb)) = p_HbO;
hemo_sim(2,ismember(info.tissue.dim.Good_Vox, perturb)) = p_HbR;
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

%% noise model
Ysim690_noise = a690*Xsim(:,1)+noise690;
Ysim850_noise = a850*Xsim(:,2)+noise850;
Xrecon690_noise = iA690*Ysim690_noise;
Xrecon850_noise = iA850*Ysim850_noise;
Xrecon_noise = [Xrecon690_noise';Xrecon850_noise'];
hemo_recon_noise = Extco^(-1)*Xrecon_noise;

%% image noise 
image_noise690 = iA690*noise690;
image_noise850 = iA850*noise850;
image_noise = Extco^(-1)*[image_noise690,image_noise850]';

end
