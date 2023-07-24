function rx_samples = simulation_channel(tx_samples, param)

h_amp = [1, 0.3]; %channel taps amplitude
h_doppler = [0, 0]; % Hz %[660,-601] % (fc*v/c)
h_delay = [0, 10]; % delay in sample index (= h_del_ns*samp_rate);
SNR = 50; %dB

%apply channel
rx_samples_temp = apply_channel_to_ofdm(tx_samples.', param.SAMP_FREQ,h_amp,h_doppler,h_delay);

% add awgn
rx_samples = awgn(rx_samples_temp,SNR);
rx_samples=rx_samples.';
end