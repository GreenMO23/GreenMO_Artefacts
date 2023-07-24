clear all;
data_dir = "C:\ag_data2\random_antenna_configs_vs_babf";
expmt_modes =["random","babf"];
num_ant = 8;

for expmt_ind=1:1:2
    expmt_mode = expmt_modes(expmt_ind);
    if(expmt_ind==1)
        test_ids = 1:120; % need more data for random configs CDF
        % Total 2^8 possibilities, the random CDF starts converging near
        % 100 trials, about for 20 trials there is some issue with LTS
        % due to random 0-1 configurations: further highlights importance
        % of the antenna selection codes
    else
        test_ids = 1:10;
        % Since BABF is deterministic, 10 trials are enough
    end
    
    [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
    avg_sinr_users_across_packets,min_sinr_users_across_packets,...
    avg_ber_users_across_packets,max_ber_users_across_packets]...
    = analyse_data_greenmo(expmt_mode,data_dir, test_ids);

    avg_evm_snrs = squeeze(mean(avg_evm_snr_users_across_packets,2));
    % avg_ber = squeeze(mean(avg_ber_users_across_tests,2))
    % discard 0 EVM SNRs, comes sometimes due to LTS correlation issues with
    % random testing
    avg_evm_snrs  = avg_evm_snrs(avg_evm_snrs ~= 0);
    % Weed out extrmely negative evm snrs: sometimes random choice is not
    % invertible, or close to not invertible
    avg_evm_snrs = avg_evm_snrs(avg_evm_snrs > -2);
    plot_handle=cdfplot(avg_evm_snrs);
    linew=3;
    plot_handle.LineWidth = linew;
    hold on

    
end

ylabel("CDF");
xlabel("SINR (dB)")
legend(["Random 0-1 antenna configs","BABF"])
