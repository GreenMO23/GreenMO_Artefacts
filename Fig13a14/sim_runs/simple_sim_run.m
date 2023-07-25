rng(1);
clear all;
meta_params.num_users = 1;
meta_params.setting_id=1;
meta_params.num_ants = 8;
rand_placement = true;
meta_params.plot_setting = true;
meta_params.bandwidth = 10e6;
meta_params.oversamp_fac = 2;
meta_params.sampling_freq = meta_params.bandwidth*meta_params.oversamp_fac;
meta_params.user1_SNR = 30; % db
meta_params.tx_pow = 20; % dbm

sim_settings_struct = load_sim_setting(meta_params, rand_placement); % load geometrical setting

channel_struct = load_channel_struct(meta_params,sim_settings_struct); % Calculate the path lengths and channel taps
 
[tx_waveforms,ofdm_tx_structs] = ofdm_tx_params(meta_params); % Create the TX waveforms

[rx_waveforms_per_user, rx_interfered_waveforms] = apply_td_chan(meta_params,tx_waveforms, channel_struct); % Convolve the channel and apply path loss to get the RX waveforms

clear channel_est_struct
% channel_mat = zeros(num_ants,num_users);

meta_params.downsamp_flag = 0;
meta_params.fig_idx = 0;
meta_params.delay_resolution = 0;
[channel_est_struct, channel_mat, orig_snr, orig_chan_pow] = synch_and_get_channels(meta_params,ofdm_tx_structs,channel_struct,rx_interfered_waveforms);
% channel_mat: num_packets*num_users*num_ants*num_subc

if(meta_params.num_users==1)
    disp("Analog beamforming")
    code_phase_shifted_rx_sigs = apply_delayed_clocks(meta_params,rx_interfered_waveforms,channel_est_struct);
    meta_params.downsamp_flag = 1;
    meta_params.fig_idx = 3;
    [channel_est_struct, channel_mat, per_ant_snr, per_ant_chan_pow] = synch_and_get_channels(meta_params,ofdm_tx_structs,channel_struct,code_phase_shifted_rx_sigs);
    summed_sigs=code_phase_shifted_rx_sigs{1};
    for ant_idx=2:1:meta_params.num_ants
        summed_sigs = summed_sigs+code_phase_shifted_rx_sigs{ant_idx};
    end
    meta_params.downsamp_flag = 1;
    meta_params.fig_idx = 5;
    [summed_snr, summed_channel_pow] = get_increased_snr(meta_params,ofdm_tx_structs,channel_struct,summed_sigs);
else
    disp("Spatial Mux.")
end

beamf_gain_chanpow = 10^((summed_channel_pow-orig_chan_pow)/20);
beamf_gain_snr_orig = summed_snr-orig_snr;
beamf_gain_snr = 10^((summed_snr-per_ant_snr)/10);

%
% eval_params = cancel_interf_equalize(ofdm_rx_struct,channel_struct,ofdm_tx_sturct.tx_bits); % Analyze the final processed RX waveforms



