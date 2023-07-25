% pre-computed data
load('./Pre_run_data/SimDataFig13a.mat')

% plot code from generate Fig13a file
lwidth=3;
figure(2)
plot_mask = [...% Algo mask description
              1,...% 1: Full Digital
              1,...% 2: FC HBF
              1,...% 3: PC HBF
              1,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              1,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x

for i=1:1:numel(plot_mask)
    if(plot_mask(i)==1)
        curr_sinrs = squeeze(sinr_mat_evals(:,i,1));
        h=cdfplot(curr_sinrs);
        h.LineWidth=lwidth;
        if(i==5 && plot_mask(5)==1) % calib trace
            h.LineStyle="--";
            h.LineWidth=lwidth-1;
        end
        hold on

    end
end


net_str = ["DBF 64RF","Fully connected","Partially connected",...
    "DBF 8 RF","N VRF N Ant (Calibration trace, should be same as DBF 1x)","GreenMO","GreenMO: FilledLN"...
    ,"DBF: 2x","DBF: 3x","DBF: 4x"];
% selected_algos = fliplr(net_str);
selected_algos = net_str(plot_mask==1);

legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("Average SINR (dB)")
ylabel("Probability")
title("")
