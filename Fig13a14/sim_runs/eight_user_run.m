function [eval_params, sinr_mat, ret_val] = eight_user_run(seed,algo_mask)

num_algos_sim = numel(algo_mask);
rng(seed);
meta_params.num_users = 4;
disable_analog_switching = false;
meta_params.setting_id=1;
meta_params.num_ants = 64;
rand_placement = 1;
meta_params.plot_setting = false;
meta_params.bandwidth = 20e6;
meta_params.oversamp_fac = meta_params.num_users;
% meta_params.oversamp_fac = 4;
meta_params.sampling_freq = meta_params.bandwidth*meta_params.oversamp_fac;
meta_params.user1_SNR = 13; % db
meta_params.tx_pow = 20; % dbm

meta_params.plot_setting=true;
sim_settings_struct = load_sim_setting(meta_params, rand_placement); % load geometrical setting
hold on
meta_params.plot_setting=false;

channel_struct = load_channel_struct(meta_params,sim_settings_struct); % Calculate the path lengths and channel taps

clear loopback tx_waveforms
[tx_waveforms, loopback,ofdm_tx_structs] = ofdm_tx_params(meta_params); % Create the TX waveforms
% loopback_noise = add_noise(loopback,meta_params);

% disp("Applying channel")

[rx_interfered_waveforms, rx_interfered_waveforms_no_noise, rx_interfered_waveforms_dig, rx_interfered_waveforms_babf] = apply_td_chan(meta_params,tx_waveforms, loopback, channel_struct); % Convolve the channel and apply path loss to get the RX waveforms

% disp("Channel applied")

clear channel_est_struct


meta_params.downsamp_flag = 0;
meta_params.fig_idx = 2;
[channel_est_struct, orig_snr,ret_val1] = synch_and_get_channels_twouser(meta_params,...
    ofdm_tx_structs,channel_struct,rx_interfered_waveforms);

if(ret_val1==-1)
    sinr_mat = zeros([num_algos_sim,3]);
    eval_params=0;
    ret_val=-1;
    return
end

[channel_est_struct_no_noise, ~,ret_val2] = synch_and_get_channels_twouser(meta_params,...
    ofdm_tx_structs,channel_struct,rx_interfered_waveforms_no_noise);
if(ret_val2==-1)
    sinr_mat = zeros([num_algos_sim,3]);
    eval_params=0;
    ret_val=-1;
    return
end


rx_interfered_waveforms_downsampled = get_downsampled_wavs(meta_params,rx_interfered_waveforms_dig);

[fully_connected_rx_sigs,partially_connected_rx_sigs] = hybrid_beamf_combiner(meta_params,...
    rx_interfered_waveforms_downsampled,channel_est_struct);

[babf_dbf1x_samples,filledln_dbf1x_samples,id_dbf1x_samples] = babf_switcher(meta_params,...
    rx_interfered_waveforms_babf,channel_est_struct,channel_est_struct_no_noise);
%     meta_params.skip_combining=true;

meta_params.skip_combining=false;

% full digital
if(algo_mask(1)==1)
    disp("Full Digital beg");
    [eval_params{1},ret_val1] = cancel_interf_equalize(meta_params,rx_interfered_waveforms_downsampled,...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    
else
    eval_params{1}=-1;
end

if(algo_mask(2)==1)
    [eval_params{2},ret_val1] = cancel_interf_equalize(meta_params,fully_connected_rx_sigs,...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("Full connected done");
else
    eval_params{2}=-1;
end

if(algo_mask(3)==1)
    [eval_params{3},ret_val1] = cancel_interf_equalize(meta_params,partially_connected_rx_sigs,...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("Partial connected done");
else
    eval_params{3}=-1;
end

if(algo_mask(4)==1)
    [eval_params{4},ret_val1] = cancel_interf_equalize(meta_params,...
        rx_interfered_waveforms_downsampled(1:meta_params.num_users),...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("User dig done");
else
    eval_params{4}=-1;
end

if(algo_mask(5)==1)
    [eval_params{5},ret_val1] = cancel_interf_equalize(meta_params,id_dbf1x_samples...
        , ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("id DBF done");
else
    eval_params{5}=-1;
end


if(algo_mask(6)==1)
    [eval_params{6},ret_val1] = cancel_interf_equalize(meta_params,babf_dbf1x_samples,...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("BABF done");
else
    eval_params{6}=-1;
end

if(algo_mask(7)==1)

    [eval_params{7},ret_val1] = cancel_interf_equalize(meta_params,filledln_dbf1x_samples,...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("Filled LN done");
else
    eval_params{7}=-1;
end

if(algo_mask(8)==1)
    [eval_params{8},ret_val1] = cancel_interf_equalize(meta_params,...
        rx_interfered_waveforms_downsampled(1:meta_params.num_users*2),...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("Dig BF 2x done");
else
    eval_params{8}=-1;
end

if(algo_mask(9)==1)
    [eval_params{9},ret_val1] = cancel_interf_equalize(meta_params,...
        rx_interfered_waveforms_downsampled(1:meta_params.num_users*3),...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("Dig BF 3x done");
else
    eval_params{9}=-1;
end


if(algo_mask(10)==1)
    [eval_params{10},ret_val1] = cancel_interf_equalize(meta_params,...
        rx_interfered_waveforms_downsampled(1:meta_params.num_users*4),...
        ofdm_tx_structs,channel_est_struct); % Analyze the final processed RX waveforms
    if(ret_val1==-1)
        sinr_mat = zeros([num_algos_sim,3]);
        eval_params=0;
        ret_val=-1;
        return
    end
    disp("Dig BF 4x done");
else
    eval_params{10}=-1;
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

    



