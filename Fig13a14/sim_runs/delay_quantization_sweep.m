rng(1);
clear all;
meta_params.num_users = 1;
meta_params.setting_id=1;
meta_params.num_ants = 8;
rand_placement = true;
meta_params.bandwidth = 10e6;
meta_params.oversamp_fac = 2;
meta_params.sampling_freq = meta_params.bandwidth*meta_params.oversamp_fac;
meta_params.user1_SNR = 20; % db
meta_params.tx_pow = 20; % dbm
meta_params.plot_setting = false;

sweep_vals = linspace(0,25,11);
sweep_nums = numel(sweep_vals);
% snr_per_ant = zeros([sweep_nums,1]);
% snr_summed = zeros([sweep_nums,1]);
beamf_gain = zeros([sweep_nums,1]);

sim_settings_struct = load_sim_setting(meta_params, rand_placement); % load geometrical setting

channel_struct = load_channel_struct(meta_params,sim_settings_struct); % Calculate the path lengths and channel taps

[tx_waveforms,ofdm_tx_structs] = ofdm_tx_params(meta_params); % Create the TX waveforms

[rx_waveforms_per_user, rx_interfered_waveforms] = apply_td_chan(meta_params,tx_waveforms, channel_struct); % Convolve the channel and apply path loss to get the RX waveforms


for sweep_idx=1:1:sweep_nums

    clear channel_est_struct
    % channel_mat = zeros(num_ants,num_users);

    meta_params.downsamp_flag = 0;
    meta_params.fig_idx = 0;
    meta_params.delay_resolution = sweep_vals(sweep_idx);
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

        beamf_gain_chanpow = 10^((summed_channel_pow-orig_chan_pow)/20);
        beamf_gain_snr_orig = summed_snr-orig_snr;
        beamf_gain_snr = 10^((summed_snr-per_ant_snr)/10);


        beamf_gain(sweep_idx) = beamf_gain_snr;
    else
        disp("Spatial Mux.")
    end

end

plot(sweep_vals,beamf_gain,'LineWidth',4)
xlabel("Clock accuracy (ns)")
ylabel("Beamforming gain (linear)")
%
% eval_params = cancel_interf_equalize(ofdm_rx_struct,channel_struct,ofdm_tx_sturct.tx_bits); % Analyze the final processed RX waveforms



