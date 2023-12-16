function [hemo_recon] = MakeFD_HemoreconfooV_neuroDOT(perturb, info, a690, a850, iA690, iA850, Extco, p_HbO, p_HbR)
%% Input the perturbation coordinates (dim1 dim2 dim3), and sensitivity
%%matrcies A and iA_psf, return the output matrix, which is the Xrecon in
%%3D space.

%Make noise perturbation
hemo_sim = zeros(2,size(a690,2));
hemo_sim(1,ismember(info.tissue.dim.Good_Vox, perturb)) = p_HbO;
hemo_sim(2,ismember(info.tissue.dim.Good_Vox, perturb)) = p_HbR;
Xsim = 1000*(Extco*hemo_sim)';

%% noise free model
Ysim690 = a690*Xsim(:,1);
Ysim850 = a850*Xsim(:,2);
Xrecon690 = iA690*Ysim690;
Xrecon850 = iA850*Ysim850;
Xrecon = [Xrecon690';Xrecon850'];
hemo_recon = Extco^(-1)*Xrecon;

end
