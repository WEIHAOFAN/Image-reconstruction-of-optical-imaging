function [noise_real690, noise_img690, noise_real850, noise_img850] = realistic_noise(r1,r2,f,rand1,rand2,a,b,c,d,g)
%% model the realistic model based on the fitted coefficients a, b, c, d, g
    if ~exist('g','var')
        a(:,1) = [0.2502, 3.933e-11, 0.6019, 1.917e-10];
        b(:,1) = [0.02913, 0.4161, 0.01052, 0.3708];
        c(:,1) = [4.625e-6,0.0105,9.685e-5,0.03573];
        d(:,1) = [0.2128,0.05585,0.1382,0.02002];
        g(:,1) = [6.769e-4,0.0013,6.785e-4,0.0013];
    end
    noise_real690 = 0.01*(a(1,1)*exp(b(1,1)*r1)+c(1,1)*exp(d(1,1)*r1))*10^(g(1,1)*(f-140)).*rand1;
    noise_img690 = (a(2,1)*exp(b(2,1)*r1)+c(2,1)*exp(d(2,1)*r1))*10^(g(2,1)*(f-140))*pi.*rand1/180;
    noise_real850 = 0.01*(a(3,1)*exp(b(3,1)*r2)+c(3,1)*exp(d(3,1)*r2))*10^(g(3,1)*(f-140)).*rand2;
    noise_img850 = (a(4,1)*exp(b(4,1)*r2)+c(4,1)*exp(d(4,1)*r2))*10^(g(4,1)*(f-140))*pi.*rand2/180;
end