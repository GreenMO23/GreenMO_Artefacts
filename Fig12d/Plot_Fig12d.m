clear all;
data_dir = "C:\ag_data2\synch_vs_unsynch";
expmt_modes = ["synch","unsynch"];
test_ids=1:1:50;

for expmt_mode = expmt_modes
    all_evm_snrs = [];
    [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
    avg_sinr_users_across_packets,min_sinr_users_across_packets,...
    avg_ber_users_across_packets,max_ber_users_across_packets,avg_evm_snr_across_users]...
    = analyse_data_greenmo(expmt_mode, data_dir, test_ids);
        
    all_evm_snrs =[all_evm_snrs avg_evm_snr_across_users];

    ph=cdfplot(10*log10(all_evm_snrs));
    ph.LineWidth = 3;
    hold on
end

ylabel("CDF");
xlabel("SINR (dB)");
legend(expmt_modes);