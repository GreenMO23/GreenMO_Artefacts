function [eval_params, sinr_mat, ret_val] = two_user_run(seed,algo_mask)

num_algos_sim = numel(algo_mask);
rng(seed);
meta_params.num_users = 2;
disable_analog_switching = false;
meta_params.setting_id=1;
meta_params.num_ants = 16;
rand_placement = true;
meta_params.plot_setting = false;
meta_params.bandwidth = 10e6;
meta_params.oversamp_fac = meta_params.num_ants;
% meta_params.oversamp_fac = 4;
meta_params.sampling_freq = meta_params.bandwidth*meta_params.oversamp_fac;
meta_params.user1_SNR = 30; % db
meta_params.tx_pow = 20; % dbm

meta_params.plot_setting=true;
sim_settings_struct = load_sim_setting(meta_params, rand_placement); % load geometrical setting
hold on
meta_params.plot_setting=false;

channel_struct = load_channel_struct(meta_params,sim_settings_struct); % Calculate the path lengths and channel taps

clear loopback tx_waveforms
[tx_waveforms, loopback,ofdm_tx_structs] = ofdm_tx_params(meta_params); % Create the TX waveforms
% loopback_noise = add_noise(loopback,meta_params);


[rx_waveforms_per_user, rx_interfered_waveforms, rx_interfered_waveforms_dig, rx_interfered_waveforms_babf] = apply_td_chan(meta_params,tx_waveforms, loopback, channel_struct); % Convolve the channel and apply path loss to get the RX waveforms
% clear test_waveform_antenna
% test_waveform_antenna{1} = resample(rx_waveforms_per_user{1,1},1,meta_params.oversamp_fac);
% test_waveform_antenna{2} = resample(rx_waveforms_per_user{2,1},1,meta_params.oversamp_fac);

clear channel_est_struct
% channel_mat = zeros(num_ants,num_users);disable_analog_switching

meta_params.downsamp_flag = 0;
meta_params.fig_idx = 0;
meta_params.delay_resolution = 0;
% channel_mat: num_packets*num_users*num_ants*num_subc


meta_params.downsamp_flag = 0;
meta_params.fig_idx = 2;
[channel_est_struct, channel_mat, orig_snr] = synch_and_get_channels_twouser(meta_params,ofdm_tx_structs,channel_struct,rx_interfered_waveforms);
% 3db insertion loss codes

rx_interfered_waveforms_downsampled = get_downsampled_wavs(meta_params,rx_interfered_waveforms_dig);

[fully_connected_rx_sigs,partially_connected_rx_sigs] = hybrid_beamf_combiner(meta_params,rx_interfered_waveforms_downsampled,channel_est_struct);

% babf4_downsampled_samples = babf4_switcher(meta_params,rx_interfered_waveforms_babf(1:meta_params.num_users*2));
switched_dbf1x_samples = babf2_switcher(meta_params,rx_interfered_waveforms_babf(1:meta_params.num_users));
clear eval_params
if(~disable_analog_switching)
    [code_phase_shifted_rx_sigs] = apply_delayed_clocks_twouser(meta_params,rx_interfered_waveforms,channel_est_struct);
    
    resampled_windowed_user_sigs = reshape_window_downsample(meta_params,code_phase_shifted_rx_sigs);
    
    [babf_dbf1x_samples,kmeans_dbf1x_samples] = babf_switcher(meta_params,rx_interfered_waveforms_babf,channel_est_struct);
    
    kmeans_2x_samples = kmeans_switcher_multirate(meta_params,rx_interfered_waveforms_babf,channel_est_struct,2);
    kmeans_4x_samples = kmeans_switcher_multirate(meta_params,rx_interfered_waveforms_dig,channel_est_struct,4);
    
    if(algo_mask(1)==1)
        meta_params.fig_idx = 4;
        [~, user_snrs_pre_combining] = observe_multiuser_channels(meta_params,ofdm_tx_structs,channel_struct,resampled_windowed_user_sigs);

        clear summed_sigs 
        for usr_idx=1:1:meta_params.num_users
            summed_sigs{usr_idx}=resampled_windowed_user_sigs{usr_idx,1};
    %         loopback{usr_idx}=tx_waveforms{usr_idx};
            for ant_idx=2:1:meta_params.num_ants
                summed_sigs{usr_idx} = summed_sigs{usr_idx} + resampled_windowed_user_sigs{usr_idx,ant_idx};
            end
        end

        meta_params.fig_idx = 6;
        meta_params.plot_setting = false;
        [digital_channel_mat, user_snrs_pre_digital, ret_val] = get_digital_channels(meta_params,ofdm_tx_structs,summed_sigs);

        if(ret_val==-1)
            sinr_mat = zeros([num_algos_sim,3]);
            eval_params=0;
            ret_val=-1;
            return
        else
            ret_val=1;
        end

        meta_params.fig_idx = 8;
        meta_params.plot_setting = false;
        meta_params.skip_combining=false;
        [eval_params{1},ret_val1] = cancel_interf_equalize(meta_params,summed_sigs, ofdm_tx_structs); % Analyze the final processed RX waveforms
        if(ret_val1==-1)
            sinr_mat = zeros([num_algos_sim,3]);
            eval_params{1}=-1;
            ret_val=-1;
            return
        end
    else
        eval_params{1}=-1;
    end
        
end
%     meta_params.skip_combining=true;

meta_params.skip_combining=false;

if(algo_mask(2)==1)
    [eval_params{2},ret_val1] = cancel_interf_equalize(meta_params,rx_interfered_waveforms_downsampled, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{2}=-1;
end

if(algo_mask(3)==1)
    [eval_params{3},ret_val1] = cancel_interf_equalize(meta_params,fully_connected_rx_sigs, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{3}=-1;
end

if(algo_mask(4)==1)
    [eval_params{4},ret_val1] = cancel_interf_equalize(meta_params,partially_connected_rx_sigs, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{4}=-1;
end

% babf4_params = cancel_interf_equalize(meta_params,babf4_downsampled_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
if(algo_mask(5)==1)
    [eval_params{5},ret_val1] = cancel_interf_equalize(meta_params,babf_dbf1x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{5}=-1;
end

if(algo_mask(6)==1)
    [eval_params{6},ret_val1] = cancel_interf_equalize(meta_params,switched_dbf1x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{6}=-1;
end    

if(algo_mask(7)==1)
    [eval_params{7},ret_val1] = cancel_interf_equalize(meta_params,rx_interfered_waveforms_downsampled(1:meta_params.num_users*3), ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{7}=-1;
end    

if(algo_mask(8)==1)
    [eval_params{8},ret_val1] = cancel_interf_equalize(meta_params,rx_interfered_waveforms_downsampled(1:meta_params.num_users*2), ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{8}=-1;
end    

if(algo_mask(9)==1)
    [eval_params{9},ret_val1] = cancel_interf_equalize(meta_params,rx_interfered_waveforms_downsampled(1:meta_params.num_users), ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{9}=-1;
end    

if(algo_mask(10)==1)
    [eval_params{10},ret_val1] = cancel_interf_equalize(meta_params,kmeans_dbf1x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{10}=-1;
end    

if(algo_mask(11)==1)
    [eval_params{11},ret_val1] = cancel_interf_equalize(meta_params,kmeans_2x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{11}=-1;
end    

if(algo_mask(12)==1)
    [eval_params{12},ret_val1] = cancel_interf_equalize(meta_params,kmeans_4x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{12}=-1;
end    



sinr_mat = zeros([num_algos_sim,3]);

for i=1:1:num_algos_sim
    if(algo_mask(i)==1)
        sinr_mat(i,1) = mean(eval_params{i}.sinr_calc,[1,2]);
        sinr_mat(i,2) = min(mean(eval_params{i}.sinr_calc,[1]));
        sinr_mat(i,3) = max(mean(eval_params{i}.sinr_calc,[1]));
    end
end

ret_val = 1;
end

%     loopback_params = cancel_interf_equalize(meta_params,test_waveform_antenna,digital_channel_mat,ofdm_tx_structs); % Analyze the final processed RX waveforms

    



