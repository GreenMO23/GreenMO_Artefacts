num_evals = 100;
% hybrid timo:1
% full digital:2
% fc hyb:3
% pc hyb:4

% babf digital 2 chains: 5
% id switch 2 chains: 6
% dig bf 3x: 7
% dig bf 2x: 8

% dig bf 1x: 9
% kmeans 1x: 10
% kmeans 2x: 11
% kmeans 4x: 12

algo_mask = [0,1,0,0,...
             1,1,1,1,...
             1,1,1,1];
sinr_mat_evals = zeros([num_evals,numel(algo_mask),3]);
clear eval_params_all
eval_idx=0;
seed_eval_idx=0;
while(eval_idx<num_evals)
    disp("Starting eval: "+num2str(eval_idx+1))
    [eval_params_all{eval_idx+1},sinr_mat_evals(eval_idx+1,:,:),ret_val] = two_user_run(seed_eval_idx,algo_mask);
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

algo_mask = [0,1,0,0,...
             0,0,0,1,...
             1,1,1,1];
         
for i=1:1:numel(algo_mask)
    if(algo_mask(i)==1)
        curr_sinrs = squeeze(sinr_mat_evals(:,i,1));
        h=cdfplot(curr_sinrs);
        h.LineWidth=lwidth;
        hold on
    end
end

h=xline(30);
h.LineWidth = 6;
net_str = ["Kmeans 8 antennas, single ADC 8*BW","Kmeans 8 antennas, single ADC, 4*BW","Kmeans 8 antennas: single ADC, 2*BW", "Digital Beamforming 2 antennas 2 RF Chains", "Digital Beamforming 4 antennas 4 RF Chains",...
"Digital Beamforming 6 RF Chains", "Identity Codes, 2 antennas 1 ADC 2*BW",...
"BABF, single ADC 2*BW", "pc hyb","fc hyb",  "Digital Beamforming 8 antennas 8 RF Chains","sw hyb"];
selected_algos = fliplr(net_str);
selected_algos = selected_algos(algo_mask==1);

legend([selected_algos,"Avg TDMA SNR"])
xlabel("Average SINR (dB)")
ylabel("CDF")
%%
figure(9)
gain_algo_mask = [0,1,0,0,...
             0,0,0,1,...
             0,1,1,1];
         
for i=1:1:numel(gain_algo_mask)
    if(gain_algo_mask(i)==1 && i~=9)
        curr_sinrs = squeeze(sinr_mat_evals(:,i,1))-squeeze(sinr_mat_evals(:,9,1));
        h=cdfplot(curr_sinrs);
        h.LineWidth=lwidth;
        hold on
    end
end

% h=xline(30);
% h.LineWidth = 6;
net_str = ["Kmeans 8 antennas, single ADC 8*BW","Kmeans 8 antennas, single ADC, 4*BW","Kmeans 8 antennas: single ADC, 2*BW", "Dig BF 2 chains", "Digital Beamforming 4 antennas 4 RF Chains",...
"Digital Beamforming 6 RF Chains", "Identity Codes, 2 antennas 1 ADC 2*BW",...
"BABF, single ADC 2*BW", "pc hyb","fc hyb",  "Digital Beamforming 8 antennas 8 RF Chains","sw hyb"];
selected_algos = fliplr(net_str);
selected_algos = selected_algos(gain_algo_mask==1);

legend([selected_algos])
xlabel("SNR Gain over 2 RF chain beamformer (dB)")
ylabel("CDF")
%%
timestamp = datestr(now,'mm_dd_HH_MM');
save("./IA_hwr_repo/simulations/Data/100evalrun_8ants_"+timestamp+".mat","sinr_mat_evals","eval_params_all","algo_mask")