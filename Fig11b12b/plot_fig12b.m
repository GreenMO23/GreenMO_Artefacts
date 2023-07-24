data_dir = "C:\ag_data2\PCB_MMSE_virt_vs_physical";
% 10 positions
posn_lists = [0,1,2,3,4,5,6,7,8,9];
% compare with mmse (dbf implemented over 4 physical chains of warp)
% pcb here refers to DBF
mode_strs = ["mmse","pcb"];

avg_mean_gput_mmse = zeros(1,3);
avg_mean_gput_osamp = zeros(1,3);
avg_std_gput_mmse = zeros(1,3);
avg_std_gput_osamp = zeros(1,3);
chan_cap = zeros(1,3);
linew = 3;

colors = ["m","b","c"];

tol = 1.5;

for num_users = 2:1:4
    tx_params_dir = data_dir+"/tx_data"+num2str(num_users)+"/";
    load(tx_params_dir+"ofdm_params_tx.mat");
    num_tx_bits = numel(ofdm_params_tx(1).tx_bits);
    tx_time(num_users-1) = (ofdm_params_tx(1).TX_NUM_SAMPS-((64)*10+32+ofdm_params_tx(1).inter_user_zeros)*num_users-ofdm_params_tx(1).N_ZERO_PAD)/10^7;
    
    sinr_cdf_list_mmse_avg=[];
    sinr_cdf_list_osamp_avg=[];
    sinr_cdf_list_mmse_min=[];
    sinr_cdf_list_osamp_min=[];
    
    evm_snr_cdf_list_mmse_avg=[];
    evm_snr_cdf_list_osamp_avg=[];
    evm_snr_cdf_list_mmse_min=[];
    evm_snr_cdf_list_osamp_min=[];
    
    ber_list_mmse_avg=[];
    ber_list_osamp_avg=[];
    ber_list_mmse_min=[];
    ber_list_osamp_min=[];
    ber_list_mmse_max=[];
    ber_list_osamp_max=[];
    capacities = [];
   

    for posn_nums = 1:1:numel(posn_lists)
        posn_index = posn_lists(posn_nums);
        expmt_dir="rx_data"+num2str(num_users)+"_"+num2str(posn_index);
        path_dir = data_dir+"/"+expmt_dir+"/";
        load(path_dir+"decode_out_pcb_mmse.mat");
    
        for trial_idx=1:1:3
            
            sinr_mat_mmse = sinr_dict{1,trial_idx};
            sinr_mat_osamp = sinr_dict{2,trial_idx};
            mean_subc_sinr_mmse = mean(sinr_mat_mmse,3);
            mean_subc_sinr_osamp = mean(sinr_mat_osamp,3);
    

            sinr_cdf_list_mmse_avg = [sinr_cdf_list_mmse_avg outlier_rej(mean(mean_subc_sinr_mmse,2).',tol)];
            sinr_cdf_list_osamp_avg = [sinr_cdf_list_osamp_avg outlier_rej(mean(mean_subc_sinr_osamp,2).',tol)];
            sinr_cdf_list_mmse_min = [sinr_cdf_list_mmse_min outlier_rej(min(mean_subc_sinr_mmse,[],2).',tol)];
            sinr_cdf_list_osamp_min = [sinr_cdf_list_osamp_min outlier_rej(min(mean_subc_sinr_osamp,[],2).',tol)];
            
            evm_snr_mat_mmse = evm_snr_dict{1,trial_idx};
            evm_snr_mat_osamp = evm_snr_dict{2,trial_idx};
            evm_snr_cdf_list_mmse_avg = [evm_snr_cdf_list_mmse_avg mean(evm_snr_mat_mmse,2).'];
            evm_snr_cdf_list_osamp_avg = [evm_snr_cdf_list_osamp_avg mean(evm_snr_mat_osamp,2).'];
            evm_snr_cdf_list_mmse_min = [evm_snr_cdf_list_mmse_min min(evm_snr_mat_mmse,[],2).'];
            evm_snr_cdf_list_osamp_min = [evm_snr_cdf_list_osamp_min min(evm_snr_mat_osamp,[],2).'];
            
            ber_mat_mmse = ber_dict{1,trial_idx};
            ber_mat_osamp = ber_dict{2,trial_idx};
            ber_list_mmse_avg = [ber_list_mmse_avg mean(ber_mat_mmse,2).'];
            ber_list_osamp_avg = [ber_list_osamp_avg mean(ber_mat_osamp,2).'];
            
            mean_tdma_snr = mean(mean(mean(tdma_snrs,4),3),1);
            capacity = sum(10^7*log10(1+mean_tdma_snr));
            capacities = [capacities capacity];
            avg_sinr;
            
        end
        clear avg_sinr avg_ber avg_snr_evm ber_dict evm_snr_dict sinr_dict
    end
    

    
    figure(1)
    plot_handle=cdfplot(10*log10(evm_snr_cdf_list_mmse_avg));
    plot_handle.LineStyle = "-.";
    plot_handle.Color = colors(num_users-1);
    plot_handle.LineWidth = linew;
    hold on
    plot_handle=cdfplot(10*log10(evm_snr_cdf_list_osamp_avg));
    plot_handle.LineStyle = "-";
    plot_handle.Color = colors(num_users-1);
    plot_handle.LineWidth = linew;
    
    mean(ber_list_mmse_avg);
    avg_mean_gput_mmse(num_users-1) = (1-mean(ber_list_mmse_avg))*num_users*num_tx_bits/tx_time(num_users-1);
    avg_mean_gput_osamp(num_users-1) = (1-mean(ber_list_osamp_avg))*num_users*num_tx_bits/tx_time(num_users-1);
    
    avg_std_gput_mmse(num_users-1) = 0.5*std(ber_list_mmse_avg)*num_users*num_tx_bits/tx_time(num_users-1);
    avg_std_gput_osamp(num_users-1) = 0.5*std(ber_list_osamp_avg)*num_users*num_tx_bits/tx_time(num_users-1);
    
    chan_cap(num_users-1) = mean(capacities);
end



figure(1)
% subplot(1,2,1)
title("SINR comparison of physical vs virtual RF chains")
ylabel("Probability")
xlabel("SINR (dB)")
ph=xline(10,'--',"Target SNR");
ph.LineWidth=3;
ph.Color="k";
xlim([-10,25])
h = zeros(5, 1);
h(1) = plot(NaN,NaN,'-m');
h(2) = plot(NaN,NaN,'-b');
h(3) = plot(NaN,NaN,'-c');
h(4) = plot(NaN,NaN,'-k');
h(5) = plot(NaN,NaN,'-.k');
legend(h, '2 Users','3 Users','4 Users','GreenMO Virtual Chains','DBF Physical Chains');


function rej_arr = outlier_rej(inp_array,tolerance)
%     [sz1,sz2] = size(inp_array);
%     std_arr = sqrt(var(inp_array,1));
%     diff_arr = abs(inp_array-repmat(std_arr,[sz1,1]));
    mean_arr = mean(inp_array);
    std_arr = sqrt(var(inp_array));
    diff_arr = abs(inp_array-mean_arr);
    rej_arr=inp_array;
    rej_arr(diff_arr>tolerance*std_arr)=[];
end