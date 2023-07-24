data_dir = "C:\ag_data2\goodput_evals";
file_name = "decode_out_fig10.mat";
posn_lists = [1,2,3,4,5,6,7,8,9,10];
avg_mean_ber_osamp = zeros(1,4);
avg_mean_gput_osamp = zeros(1,4);
avg_std_gput_osamp = zeros(1,4);
chan_cap = zeros(1,4);
tdma_snrs_all = zeros(1,4);
linew = 3;

colors = ["r","m","b","c"];
num_ants_list = [4,6,8];

num_users = 4;
tol = 1.5;
tdma_snr_est_evm = 16;
target_3x_evm=10*log10(2^((3/4)*log2(1+10^(tdma_snr_est_evm/10)))-1);
target_2x_evm=10*log10(2^((2/4)*log2(1+10^(tdma_snr_est_evm/10)))-1);

for num_ant_idx = 1:1:numel(num_ants_list)
    
    tx_params_dir = data_dir+"/tx_data"+num2str(num_users)+"/";
    load(tx_params_dir+"ofdm_params_tx.mat");
    tx_time = (ofdm_params_tx(1).TX_NUM_SAMPS-((64)*8+ofdm_params_tx(1).inter_user_zeros)*num_users-ofdm_params_tx(1).N_ZERO_PAD)/10^7;
    num_tx_bits = numel(ofdm_params_tx(1).tx_bits);
    
    
    sinr_cdf_list_osamp_avg=[];
    sinr_cdf_list_osamp_min=[];
    
    evm_snr_cdf_list_osamp_avg=[];
    evm_snr_cdf_list_osamp_min=[];
    
    ber_list_osamp_avg=[];
    ber_list_osamp_min=[];
    ber_list_osamp_max=[];
    
    gput_cdf_list = [];
    
    capacities = [];
    tdma_snrs_list = [];
    num_ant = num_ants_list(num_ant_idx);
    

    for posn_nums = 1:1:numel(posn_lists)
        posn_index = posn_lists(posn_nums);

        expmt_str="rx_data"+num2str(num_ant)+"_"+num2str(posn_index);
        path_dir = data_dir+"/"+expmt_str+"/";
        
        load(path_dir+file_name);
        
        for trial_idx=1:1:3
            
            sinr_mat_osamp = sinr_dict{1,trial_idx};
            mean_subc_sinr_osamp = mean(sinr_mat_osamp,3);
            avg_snrs_across_users = mean(mean_subc_sinr_osamp,2).';
            min_snr_across_users = min(mean_subc_sinr_osamp,[],2).';
            % outlier rejection
            std_snr = std(avg_snrs_across_users);
            mean_snr = mean(avg_snrs_across_users);
            avg_snrs_across_users(abs(avg_snrs_across_users-mean_snr)>tol*std_snr)=[];
            
            sinr_cdf_list_osamp_avg = [sinr_cdf_list_osamp_avg avg_snrs_across_users];
            sinr_cdf_list_osamp_min = [sinr_cdf_list_osamp_min min_snr_across_users];
            
            evm_snr_mat_osamp = evm_snr_dict{1,trial_idx};
            avg_snrs_across_users = mean(evm_snr_mat_osamp,2).';
            min_snr_across_users = min(evm_snr_mat_osamp,[],2).';
            % outlier rejection
            std_snr = std(avg_snrs_across_users);
            mean_snr = mean(avg_snrs_across_users);
            avg_snrs_across_users(abs(avg_snrs_across_users-mean_snr)>tol*std_snr)=[];
            evm_snr_cdf_list_osamp_avg = [evm_snr_cdf_list_osamp_avg avg_snrs_across_users];
            evm_snr_cdf_list_osamp_min = [evm_snr_cdf_list_osamp_min min_snr_across_users];
            
            ber_mat_osamp = ber_dict{1,trial_idx};
            avg_ber_across_users = mean(ber_mat_osamp,2).';
            std_ber = std(avg_ber_across_users);
            mean_ber = mean(avg_ber_across_users);
            avg_ber_across_users(abs(avg_ber_across_users-mean_ber)>tol*std_ber)=[];
            ber_list_osamp_avg = [ber_list_osamp_avg avg_ber_across_users];
            
            curr_gput = (1-avg_ber_across_users)*48;
            gput_cdf_list = [gput_cdf_list curr_gput];
            
            
            mean_tdma_snr = mean(tdma_snrs,[1,3,4]);
            capacity = sum(10^7*log10(1+mean_tdma_snr));
            capacities = [capacities capacity];
            tdma_snrs_list = [tdma_snrs_list mean(mean_tdma_snr)];
            
        end
        clear avg_sinr avg_ber avg_snr_evm ber_dict evm_snr_dict sinr_dict
    end
    
    figure(1)
    plot_handle=cdfplot(10*log10(sinr_cdf_list_osamp_avg));
    plot_handle.LineStyle = "-";
    plot_handle.Color = colors(num_ant_idx);
    plot_handle.LineWidth = linew;
    hold on
    
    figure(2)
    plot_handle=cdfplot(10*log10(evm_snr_cdf_list_osamp_avg));
    plot_handle.LineStyle = "-";
    plot_handle.Color = colors(num_ant_idx);
    plot_handle.LineWidth = linew;
    hold on
    
    figure(3)
    plot_handle=cdfplot(gput_cdf_list);
    mean(gput_cdf_list)
    std(gput_cdf_list)
    plot_handle.LineStyle = "-";
    plot_handle.Color = colors(num_ant_idx);
    plot_handle.LineWidth = linew;
    hold on
    
    avg_mean_gput_osamp(num_ant_idx) = (1-mean(ber_list_osamp_avg))*num_users*num_tx_bits/tx_time;
    
    avg_std_gput_osamp(num_ant_idx) = 0.5*std(ber_list_osamp_avg)*num_users*num_tx_bits/tx_time;
    
    chan_cap(num_ant_idx) = mean(capacities);
    
    tdma_snrs_all(num_ant_idx) = 10*log10(mean(tdma_snrs_list));
    
    avg_mean_ber_osamp(num_ant_idx) = mean(ber_list_osamp_avg);
end


legend_str = ["4 antenna GreenMO","6 antenna GreenMO","8 antenna GreenMO"];
figure(1)
title("LTS SINR (mean) after cancellation CDFs (approx +5 offset from EVM method)")
ylabel("Probability")
xlabel("SINR (dB)")
ph=xline(15,'--',"Target SNR");
ph.LineWidth=3;
ph.Color="k";
xlim([0,35])
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-r');
h(2) = plot(NaN,NaN,'-m');
h(3) = plot(NaN,NaN,'-b');
legend(h, legend_str);

figure(2)
title("EVM SNR (mean) after cancellation CDFs")
ylabel("Probability")
xlabel("SINR (dB)")
ph=xline(15,'--',"Target SINR");
ph.LineWidth=3;
ph.Color="k";
ph.LabelVerticalAlignment = 'bottom';
xlim([-10,20])
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-r');
h(2) = plot(NaN,NaN,'-m');
h(3) = plot(NaN,NaN,'-b');
legend(h, legend_str);


figure(3)
title("Goodput CDFs")
ylabel("Probability")
xlabel("Goodput (Mbps)")
ph=xline(48,'--',"Max Goodput: 48 Mbps");
ph.LineWidth=5;
ph.Color="k";
xlim([10,60])
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-r');
h(2) = plot(NaN,NaN,'-m');
h(3) = plot(NaN,NaN,'-b');
legend(h, legend_str);
