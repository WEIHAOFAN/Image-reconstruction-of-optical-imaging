mc_time2 = 0:1e-10:5e-9-1e-11;
det_elem = mesh.elements(mesh.source.int_func(114,1),:);
mc_tpsf2 = squeeze(flux.data(det_elem(1),:));
figure
semilogy(mc_time2,abs(mc_tpsf2(1,:)),'o-','LineWidth',1.5)