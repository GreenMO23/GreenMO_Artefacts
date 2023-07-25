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

user_sweeps = 2:2:8;
sinr_mat_evals = zeros([num_evals,numel(user_sweeps),numel(algo_mask),4]);
clear eval_params_all
eval_idx=0;
seed_eval_idx=1;
while(eval_idx<num_evals)
    disp("Starting eval: "+num2str(eval_idx+1))
    [eval_params_all{eval_idx+1},sinr_mat_evals(eval_idx+1,:,:,:),ret_val] = user_sweep_run(seed_eval_idx,algo_mask);
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


% rfc_pow = 1.2*1000;
% var_pow_const = (((1890)/64)); %4200*0.45
% cool_pow = 0.05*4200;

% 45,30,20,5
% rfc_pow = 1.2*1000;
bb_pow_const = (((1260)/64)); %4200*0.3
rf_pow_const = (((840)/64)); %4200*0.2
% aa_pow_const_pc = 
cool_pow = 0.05*4200;


fixed_pow = [ ...% Algo mask description
              cool_pow+1260+840,...% 1: Full Digital
              cool_pow+rf_pow_const*8*8,...% 2: FC HBF
              cool_pow+rf_pow_const*8,...% 3: PC HBF
              cool_pow,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              cool_pow+rf_pow_const,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x

var_pow = [ ...% Algo mask description
              0,...% 1: Full Digital
              bb_pow_const,...% 2: FC HBF
              bb_pow_const,...% 3: PC HBF
              rf_pow_const+bb_pow_const,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              bb_pow_const,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x



for i=1:1:numel(plot_mask)
    if(plot_mask(i)==1)
        % mean sinrs: evals*sweeps
        
        curr_mean_sinrs_across_users = sinr_mat_evals(:,:,i,1);
        curr_mean_std_across_users = sinr_mat_evals(:,:,i,4);
%         curr_spread_sinrs_across_users = sinr_mat_evals(:,:,i,3)-sinr_mat_evals(:,:,i,2);
        mean_sinr_across_evals = squeeze(mean(curr_mean_sinrs_across_users,1));
%         mean_spread = squeeze(mean(curr_spread_sinrs_across_users,1));
        mean_std_across_evals = squeeze(mean(curr_mean_std_across_users,1));
        figure(1)
        % h=plot(user_sweeps, mean_sinr_across_evals);
        h=errorbar(user_sweeps,mean_sinr_across_evals,mean_std_across_evals/2);
        
        h.LineWidth=lwidth;
%         h.LineStyle="-.";

%         if(i==6 && plot_mask(5)==1) % calib trace
% %             h.LineStyle="--";
% %             h.LineWidth=lwidth-1;
%         end
        hold on

        figure(2)
        curr_capacity = 100*log2(1+10.^(mean_sinr_across_evals/10)).*user_sweeps;
        net_pow = fixed_pow(i)+var_pow(i).*user_sweeps;
        h = plot(user_sweeps,curr_capacity./net_pow);
        h.LineStyle="--";
        h.Marker="*";
        h.MarkerSize=10;
        h.LineWidth=lwidth;
        hold on

        figure(3)
        curr_capacity = 100*log2(1+10.^(mean_sinr_across_evals/10)).*user_sweeps;
%         net_pow = fixed_pow(i)+var_pow(i).*user_sweeps;
        h = plot(user_sweeps,curr_capacity/100);
        h.LineStyle="--";
        h.Marker="*";
        h.MarkerSize=10;
        h.LineWidth=lwidth;
        hold on
        
    end
end


net_str = ["DBF 64RF","Fully connected","Partially connected",...
    "DBF 8 RF","N VRF N Ant (Calib, should be same as DBF 1x)","GreenMO","GreenMO: FilledLN"...
    ,"DBF: 2x","DBF: 3x","DBF: 4x"];
% selected_algos = fliplr(net_str);
selected_algos = net_str(plot_mask==1);

figure(1)
legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("Number of users")
ylabel("SINR (dB)")
title("")
grid on

figure(2)
legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("Number of users")
ylabel("Energy Efficiency (Mb per Joule)")
title("")
grid on

figure(3)
legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("Number of users")
ylabel("Spectral Efficiency (Bits per Hz)")
title("")
grid on
timestamp = datestr(now,'mm_dd_HH_MM');
%%
save("./IA_hwr_repo/simulations/Data/10_unsynch_evalrun_64ants_"+timestamp+".mat","sinr_mat_evals","eval_params_all","algo_mask")