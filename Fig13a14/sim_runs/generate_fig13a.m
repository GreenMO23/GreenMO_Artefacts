num_evals = 100;

algo_mask = [...% Algo mask description
              1,...% 1: Full Digital
              1,...% 2: FC HBF
              1,...% 3: PC HBF
              1,...% 4: UP Dig
              1,...% 5: id DBF (Calib)
              1,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];% 10: DBF 4x

sinr_mat_evals = zeros([num_evals,numel(algo_mask),3]);
clear eval_params_all
eval_idx=0;
seed_eval_idx=0;
while(eval_idx<num_evals)
    disp("Starting eval: "+num2str(eval_idx+1))
    [eval_params_all{eval_idx+1},sinr_mat_evals(eval_idx+1,:,:),ret_val] = eight_user_run(seed_eval_idx,algo_mask);
    if (ret_val==-1)
        seed_eval_idx=seed_eval_idx+1;
        continue
    else
        seed_eval_idx=seed_eval_idx+1;
        eval_idx=eval_idx+1;
    end
end

%%
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

timestamp = datestr(now,'mm_dd_HH_MM');
%%
save("./IA_hwr_repo/simulations/Data/10_unsynch_evalrun_64ants_"+timestamp+".mat","sinr_mat_evals","eval_params_all","algo_mask")