function syms_eq_mat = rx_equalize_signal(rx_vec_cfo_corr,rx_H_est,param)
% Take entire rx signal, get location of one packet using lts_ind.
% Then take one packet and equalize the signal with channel rx_H_est
% Return the equalized packet syms_eq_mat of size 64x100 for 64subs and
% 100symbols



payload_ind = 1+((param.NUM_LTS+0.5)*param.N_SC+param.inter_user_zeros)*(param.num_users);

% number of original samples 8000 = 80*100
tx_num_origsamp = (param.N_SC+param.CP_LEN)*param.N_OFDM_SYMS;

%packet index range: rx_lts+160 to rx_lts+160+8000-1
rx_pkt_range = payload_ind:payload_ind+tx_num_origsamp-1;

% Extract the received packet (remove preamble as well)
payload_vec = rx_vec_cfo_corr(rx_pkt_range);% lts_ind_n+160)

% Reshape
payload_mat = reshape(payload_vec, (param.N_SC+param.CP_LEN), param.N_OFDM_SYMS);

% Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
payload_mat_noCP = payload_mat(param.CP_LEN-param.FFT_OFFSET+[1:param.N_SC], :);

% Take the FFT
syms_f_mat = fft(payload_mat_noCP, param.N_SC, 1);

% Equalize (zero-forcing, just divide by complex chan estimates)

rx_H_est=rx_H_est(:);
syms_eq_mat = syms_f_mat ./ (repmat(rx_H_est, 1, param.N_OFDM_SYMS)+1e-6);

end