function [eval_params, sinr_mat, ret_val] = two_user_ant_sweep(seed,algo_mask, num_ants)

num_algos_sim = numel(algo_mask);
rng(seed);
meta_params.num_users = 2;
disable_analog_switching = false;
meta_params.setting_id=1;
meta_params.num_ants = num_ants;
rand_placement = true;
meta_params.plot_setting = false;
meta_params.bandwidth = 10e6;
% meta_params.oversamp_fac = meta_params.num_ants;
meta_params.oversamp_fac = meta_params.num_users;
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
[channel_est_struct, channel_mat, orig_snr, ret_val1] = synch_and_get_channels_twouser(meta_params,ofdm_tx_structs,channel_struct,rx_interfered_waveforms_dig);
% 3db insertion loss codes
if(ret_val1==-1)
    sinr_mat = zeros([num_algos_sim,3]);
    eval_params=0;
    ret_val=-1;
    return
end
    
rx_interfered_waveforms_downsampled = get_downsampled_wavs(meta_params,rx_interfered_waveforms_dig);

[fully_connected_rx_sigs,partially_connected_rx_sigs] = hybrid_beamf_combiner(meta_params,rx_interfered_waveforms_downsampled,channel_est_struct);

% babf4_downsampled_samples = babf4_switcher(meta_params,rx_interfered_waveforms_babf(1:meta_params.num_users*2));
switched_dbf1x_samples = babf2_switcher(meta_params,rx_interfered_waveforms_dig(1:meta_params.num_users));
clear eval_params


[babf_dbf1x_samples,kmeans_dbf1x_samples] = babf_switcher(meta_params,rx_interfered_waveforms_dig,channel_est_struct);
    
meta_params.skip_combining=false;

% Full Digital, all antennas
if(algo_mask(1)==1)
    [eval_params{1},ret_val1] = cancel_interf_equalize(meta_params,rx_interfered_waveforms_downsampled, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{1}=-1;
end

% Fully connected PSs
if(algo_mask(2)==1)
    [eval_params{2},ret_val1] = cancel_interf_equalize(meta_params,fully_connected_rx_sigs, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{2}=-1;
end

% Partially connected PSs
if(algo_mask(3)==1)
    [eval_params{3},ret_val1] = cancel_interf_equalize(meta_params,partially_connected_rx_sigs, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{3}=-1;
end

% All antennas + BABF
if(algo_mask(4)==1)
    [eval_params{4},ret_val1] = cancel_interf_equalize(meta_params,babf_dbf1x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{4}=-1;
end

% All antennas + Id switching
if(algo_mask(5)==1)
    [eval_params{5},ret_val1] = cancel_interf_equalize(meta_params,switched_dbf1x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{5}=-1;
end    

% All antennas + lmeans switching
if(algo_mask(6)==1)
    [eval_params{6},ret_val1] = cancel_interf_equalize(meta_params,kmeans_dbf1x_samples, ofdm_tx_structs); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
else
    eval_params{6}=-1;
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

    



