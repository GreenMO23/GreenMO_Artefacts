%July 8: Copied from OTFS_trx_v2.m
function plot_constellation(M, tx_syms_mat,rx_sim_plot_cfo)

plot(rx_sim_plot_cfo,'ro','MarkerSize',1);
hold on;grid on;
if(M==2)
    plot(tx_syms_mat(:),zeros(1,length(tx_syms_mat(:))),'bo');
else
    plot(tx_syms_mat(:),'bo');
end
axis square; axis(1.5*[-1 1 -1 1]);
title(['Tx and Rx Constellations: M=',num2str(M)])