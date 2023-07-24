function [rx_H_est,rx_vec_cfo_corr,params]= extractChannel_singuser(raw_rx_dec,params,user_idx)
% Take rx signal and lts ind given by pack_detect() function
% Set few variables in param that are used in this function.
% Finally extract the LTS and get channel estimation as CSI
% also returns CFO corrected waveform (same size as raw_rx_dec)
% CFO value is returned in param.rx_cfo_est_lts

lts_ind = 1+((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)*(user_idx);
lts_ind_range = (2-1+0.5)*params.N_SC+1:(2+0.5)*params.N_SC; % 97:160 for N_SC=64

%Extract LTS (not yet CFO corrected)
rx_lts = raw_rx_dec(lts_ind : lts_ind+(2+0.5)*params.N_SC-1); %leaving CFO detec undisturbed with numLTS
rx_lts1 = rx_lts(-params.N_SC+-params.FFT_OFFSET + lts_ind_range);
rx_lts2 = rx_lts(-params.FFT_OFFSET + lts_ind_range);

%Calculate coarse CFO est
lts_xcor = xcorr(rx_lts2,rx_lts1);
rx_cfo_est_lts = angle(lts_xcor(params.N_SC));
rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*params.N_SC);

% CFOestimate = rx_cfo_est_lts*params.SAMP_FREQ;


params.rx_cfo_est_lts = rx_cfo_est_lts;
% Apply CFO correction to raw Rx waveform
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(raw_rx_dec)-1]);
rx_vec_cfo_corr = raw_rx_dec(:) .* rx_cfo_corr_t(:);

% Reextract LTS of CFO corrected signal.
% rx_lts1 = rx_lts(-params.N_SC+-params.FFT_OFFSET + lts_ind_range);
% rx_lts2 = rx_lts(-params.FFT_OFFSET + lts_ind_range);
lts_len = (params.NUM_LTS+0.5)*params.N_SC;
% rx_lts = zeros(params.num_users,lts_len);

rx_H_est_all_lts = zeros(params.NUM_LTS,params.N_SC);

rx_lts = rx_vec_cfo_corr(lts_ind : lts_ind+lts_len);
lts_ind_range = (params.NUM_LTS-1+0.5)*params.N_SC+1:(params.NUM_LTS+0.5)*params.N_SC; % 97:160 for N_SC=64
    
for i = 1:1:params.NUM_LTS
    rx_lts_indiv = rx_lts(-params.N_SC*(i-1)-params.FFT_OFFSET + lts_ind_range);
%         xline(lts_ind-params.N_SC*(i-1)-params.FFT_OFFSET+max(lts_ind_range),'r--')
    rx_lts_f = fft(rx_lts_indiv);
    rx_H_est_all_lts(i,:) = params.lts_f(:) .* rx_lts_f(:);
end

rx_H_est = squeeze(mean(rx_H_est_all_lts,1));
    
end