% posn_index = 0;
% folder_name = "8ant_multipos_QAM16_34r"
% test_ids=1:1:3;
% avg_evm_snr_users_across_packets = zeros([num_tests,num_users]);
% min_evm_snr_users_across_packets = zeros([num_tests,num_users]);
% avg_sinr_users_across_packets = zeros([num_tests,num_users]);
% min_sinr_users_across_packets = zeros([num_tests,num_users]);
% avg_ber_users_across_packets = zeros([num_tests,num_users]);
% max_ber_users_across_packets = zeros([num_tests,num_users]);
function [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
           avg_sinr_users_across_packets,min_sinr_users_across_packets,...
           avg_ber_users_across_packets,max_ber_users_across_packets]...
           = analyse_data_fdma(posn_index,folder_name, test_ids)


num_ant = 1;
root_dir = "/media/agrim/57810C9C5890F0C7/";
savpath_dir = root_dir+"conf_room_fdma/"+folder_name;
num_users = 1;
tx_params_dir = savpath_dir+"/tx_data"+num2str(num_users)+"/";
load(tx_params_dir+"ofdm_params_tx.mat");
expmt_str="rx_data"+num2str(num_ant)+"_"+num2str(posn_index);
path_dir = savpath_dir+"/"+expmt_str+"/";

num_bands = 4;

%% Parameters and Flags
clear op_struct_fdma rx_samps raw_rx_samps rx_samples_resamped ofdm_params_out 
clear avg_sinr avg_ber avg_snr_evm ber_dict evm_snr_dict sinr_dict
                    
ofdm_params = wcsng_ofdm_param_gen(64);
% For OFDM RX
count = 2e6; % roughly 40e6 is 1sec data for 40.96MSPS receiver.
start = 0.5e6; % remove initial 1e6 samples because they might be corrupted
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 16;
ofdm_params.NUM_LTS = 10;
num_lts_to_use = 10;
ofdm_params.inter_user_zeros = 512;
ofdm_params.num_users = num_users;
ofdm_params.enable_channel_codes = true;


[tx_samples,ofdm_params] = ofdm_tx(ofdm_params); % Take tx packet and params from OFDM_tx.m
plot_debug=false;
disp_plots=false;
samp_fac = 4;

linew=3;

min_evm_snrs=[];
num_tests = numel(test_ids);
% mean_lts_snr = zeros(num_bands,numel(test_ids));
% mean_evm_snr = zeros(num_bands,numel(test_ids));
avg_evm_snr_users_across_packets = zeros([num_tests,num_bands]);
min_evm_snr_users_across_packets = zeros([num_tests,num_bands]);
avg_sinr_users_across_packets = zeros([num_tests,num_bands]);
min_sinr_users_across_packets = zeros([num_tests,num_bands]);
avg_ber_users_across_packets = zeros([num_tests,num_bands]);
max_ber_users_across_packets = zeros([num_tests,num_bands]);

for band_idx=1:1:num_bands
    disp("Test idx: Band: "+num2str(band_idx));
    for test_index = 1:1:num_tests
        % Load data
        test_idx = test_ids(test_index);
%         disp("Test idx: "+num2str(test_idx)+", Users: "+num2str(band_idx));
        fname = path_dir+"rx_samps_fdma_"+num2str(test_idx);
        load(fname);


        %% Downsample and phase align appropriately
        
        rx_samples_resamped = resample(squeeze(raw_rx_samps_fdma(band_idx,:)),1,samp_fac);    
        
        %% decode
        ofdm_params.DO_DECODE = 1;
        ofdm_params.DO_APPLY_CFO_CORRECTION = 1;
        ofdm_params.DO_APPLY_SFO_CORRECTION = 1;
        ofdm_params.DO_APPLY_PHASE_ERR_CORRECTION = 1;

        [op_struct_fdma{band_idx,test_idx}, ofdm_params_out] = ofdm_rx_siso(rx_samples_resamped, ofdm_params,ofdm_params_tx);
        if(isempty(op_struct_fdma{band_idx,test_idx}))
            avg_evm_snr_users_across_packets(test_index,band_idx) = -1;
            min_evm_snr_users_across_packets(test_index,band_idx) = -1;
            avg_sinr_users_across_packets(test_index,band_idx) = -1;
            min_sinr_users_across_packets(test_index,band_idx) = -1;
            avg_ber_users_across_packets(test_index,band_idx) = -1;
            max_ber_users_across_packets(test_index,band_idx) = -1;
            ber_dict_fdma{test_index,band_idx} = -1;
            evm_snr_dict_fdma{test_index,band_idx} = -1;
            sinr_dict_fdma{test_index,band_idx} = -1;
            disp("LTS Corr not found, exiting with zeros")
            continue
        end
        avg_evm_snr_users_across_packets(test_index,band_idx) = mean(op_struct_fdma{band_idx,test_idx}.snr);
        min_evm_snr_users_across_packets(test_index,band_idx) = min(op_struct_fdma{band_idx,test_idx}.snr);
        avg_sinr_users_across_packets(test_index,band_idx) = mean(op_struct_fdma{band_idx,test_idx}.meansnr);
        min_sinr_users_across_packets(test_index,band_idx) = min(op_struct_fdma{band_idx,test_idx}.meansnr);
        avg_ber_users_across_packets(test_index,band_idx) = mean(op_struct_fdma{band_idx,test_idx}.berratio);
        max_ber_users_across_packets(test_index,band_idx) = max(op_struct_fdma{band_idx,test_idx}.berratio);
        
        ber_dict_fdma{test_index,band_idx} = op_struct_fdma{band_idx,test_idx}.berratio;
        evm_snr_dict_fdma{test_index,band_idx} = op_struct_fdma{band_idx,test_idx}.snr;
        sinr_dict_fdma{test_index,band_idx} = op_struct_fdma{band_idx,test_idx}.meansnr;
        
    end
end

save(path_dir+"/decode_out_fdma_.mat",'ber_dict_fdma','evm_snr_dict_fdma','sinr_dict_fdma');

end

