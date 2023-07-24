function [rx_H_est,rx_vec_cfo_corr,params]= extractChannel_siso(raw_rx_dec,lts_ind,params)
% Take rx signal and lts ind given by pack_detect() function
% Set few variables in param that are used in this function.
% Finally extract the LTS and get channel estimation as CSI
% also returns CFO corrected waveform (same size as raw_rx_dec)
% CFO value is returned in param.rx_cfo_est_lts

lts_ind_range = (2-1+0.5)*params.N_SC+1:(2+0.5)*params.N_SC; % 97:160 for N_SC=64
fft_offset = 0;
if(params.DO_APPLY_CFO_CORRECTION)
    %Extract LTS (not yet CFO corrected)
    rx_lts = raw_rx_dec(lts_ind : lts_ind+(2+0.5)*params.N_SC-1); %leaving CFO detec undisturbed with numLTS
    rx_lts1 = rx_lts(-params.N_SC+-fft_offset+ lts_ind_range);
    rx_lts2 = rx_lts(-fft_offset + lts_ind_range);
    
    %Calculate coarse CFO est
    lts_xcor = xcorr(rx_lts2,rx_lts1);
    rx_cfo_est_lts = angle(lts_xcor(params.N_SC));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*params.N_SC);
    %Q:Ish-why dividing by 64 here? because LTS are N_SC apart.
    CFOestimate = rx_cfo_est_lts*params.SAMP_FREQ;
else
    rx_cfo_est_lts = 0;
end
params.rx_cfo_est_lts = rx_cfo_est_lts;
% Apply CFO correction to raw Rx waveform
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(raw_rx_dec)-1]);
rx_vec_cfo_corr = raw_rx_dec(:) .* rx_cfo_corr_t(:);

% Reextract LTS of CFO corrected signal.
% rx_lts1 = rx_lts(-params.N_SC+-params.FFT_OFFSET + lts_ind_range);
% rx_lts2 = rx_lts(-params.FFT_OFFSET + lts_ind_range);
rx_lts = rx_vec_cfo_corr(lts_ind : lts_ind+(params.NUM_LTS+0.5)*params.N_SC-1); 
lts_ind_range = (params.NUM_LTS-1+0.5)*params.N_SC+1:(params.NUM_LTS+0.5)*params.N_SC; % 97:160 for N_SC=64
rx_H_est_f = zeros(params.NUM_LTS,params.N_SC);
for i = 1:1:params.NUM_LTS
    rx_lts_indiv = rx_lts(-params.N_SC*(i-1)-params.FFT_OFFSET + lts_ind_range);
    rx_lts_f = fft(rx_lts_indiv);
    rx_H_est_f(i,:) = params.lts_f(:) .* rx_lts_f(:);
end
    
    
%take fft
% rx_lts1_f = fft(rx_lts1);
% rx_lts2_f = fft(rx_lts2);

% Calculate channel estimate from average of 2 training symbols
% rx_H_est = params.lts_f(:) .* (rx_lts1_f(:) + rx_lts2_f(:))/2;
rx_H_est = mean(rx_H_est_f(9:10,:),1);
% rx_H_est = mean(rx_H_est_f,1);
%SNR estimate: this is better than EVM method and gives ~5db higher SNR
%estimate. 
% meanvalue = (rx_lts1_f(:) + rx_lts2_f(:))/2;
% varvalue = var([rx_lts1_f(:),rx_lts2_f(:)],0,2);

meanvalue = mean(abs(rx_H_est_f),1);
varvalue = var(abs(rx_H_est_f),1);

snrestimate = meanvalue./sqrt(varvalue + 1e-20);
meansnr = db(mean(snrestimate(params.SC_IND_DATA)));
params.meansnr= meansnr;

if(params.plot_flag)
    figure(params.cf);clf;
    plot(db(snrestimate),'linewidth',2)
    hold on; yline(meansnr,'r','linewidth',3);
    % plot([1 length(evm_mat(:))], 100*[aevms(1), aevms(1)],'r','LineWidth',4)
    myAxis = axis;
    h = text(round(.05*params.N_SC), meansnr+ .1*(myAxis(4)-myAxis(3)), sprintf('Effective SNR: %.1f dB', meansnr));
    set(h,'Color',[1 0 0])
    set(h,'FontWeight','bold')
    set(h,'FontSize',10)
    set(h,'EdgeColor',[1 0 0])
    set(h,'BackgroundColor',[1 1 1])
    xlabel('subcarrier index'); ylabel('SNR (dB)');
end
end