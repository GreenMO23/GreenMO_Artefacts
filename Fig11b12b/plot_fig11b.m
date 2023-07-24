data_dir = "C:\ag_data2\PCB_MMSE_virt_vs_physical";

posn_lists = [0,1,2,3,4,5,6,7,8,9];
mode_strs = ["pcb"];

avg_mean_gput_mmse = zeros(1,3);
avg_mean_gput_osamp = zeros(1,3);
avg_std_gput_mmse = zeros(1,3);
avg_std_gput_osamp = zeros(1,3);
chan_cap = zeros(1,3);
linew = 3;

colors = ["r","m","b","c"];

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
            
        end
        clear avg_sinr avg_ber avg_snr_evm ber_dict evm_snr_dict sinr_dict
    end
    
    figure(1)
    goodput = (1-ber_list_osamp_avg)*num_users*num_tx_bits/tx_time(num_users-1);
    plot_handle=cdfplot(goodput/10^6);
    plot_handle.LineStyle = "-";
    plot_handle.Color = colors(num_users-1);
    plot_handle.LineWidth = linew;

    hold on
    
    mean(ber_list_mmse_avg);
    avg_mean_gput_mmse(num_users-1) = (1-mean(ber_list_mmse_avg))*num_users*num_tx_bits/tx_time(num_users-1);
    avg_mean_gput_osamp(num_users-1) = (1-mean(ber_list_osamp_avg))*num_users*num_tx_bits/tx_time(num_users-1);
    
    avg_std_gput_mmse(num_users-1) = 0.5*std(ber_list_mmse_avg)*num_users*num_tx_bits/tx_time(num_users-1);
    avg_std_gput_osamp(num_users-1) = 0.5*std(ber_list_osamp_avg)*num_users*num_tx_bits/tx_time(num_users-1);
    
    chan_cap(num_users-1) = mean(capacities);
end

figure(1)
% Single Stream Goodput
goodput = 12; 
plot_handle=cdfplot(goodput);
plot_handle.LineStyle = "-";
plot_handle.Color = colors(4);
plot_handle.LineWidth = linew;
hold on


figure(1)
title("Multi-stream goodput CDFs")
ylabel("Probability")
xlabel("Goodput (Mbps)")
ph=xline(12,'--',"1 Stream max goodput: 12 Mbps");
ph=xline(24,'--',"2 Streams max goodput: 24 Mbps");
ph=xline(36,'--',"3 Streams max goodput: 36 Mbps");
ph=xline(48,'--',"4 Streams max goodput: 48 Mbps");
xlim([0,50])
% ph.LineWidth=3;
% ph.Color="k";
% xlim([0,40])
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'-r');
h(2) = plot(NaN,NaN,'-m');
h(3) = plot(NaN,NaN,'-b');
h(4) = plot(NaN,NaN,'-c');
% h(6) = plot(NaN,NaN,'-.k');

legend(h,'1 Stream', '2 Streams','3 Streams','4 Streams');


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