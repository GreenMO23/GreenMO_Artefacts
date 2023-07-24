function [op_struct, params] = ofdm_rx(rx_vec_air, params,first_lts_ind...
                                    ,lts_thresh)
%OFDM_RX OFDM SISO receiver, from WARPLab, Generalized to to any number of subcarriers


if nargin == 4
    lts_predetec = false;
elseif nargin == 3
    lts_predetec = true;
elseif nargin == 2
    lts_predetec = false;
    lts_thresh = 0.6;
else
    disp("Wrong param usage")
    op_stuct =[];
    params=[];
    return
end

raw_rx_dec = rx_vec_air;
% Correlate for LTS
if(~lts_predetec)
    [first_lts_ind] = pack_detect(raw_rx_dec,params,lts_thresh);
    op_struct.lts_ind = first_lts_ind;
else
    op_struct.lts_ind = first_lts_ind;
end
    

%% NEW THINGS

if(first_lts_ind<0)
    op_struct= [];
    params = [];
    return
end

% length_rx_samples = length(raw_rx_dec)-first_lts_ind-params.TX_NUM_SAMPS+1;%-2*1.5*params.TX_NUM_SAMPS; %max index for the start of a packet
first_lts_ind = op_struct.lts_ind(1,1)-((params.NUM_LTS)*params.N_SC);
op_struct.user_lts_offsets = op_struct.lts_ind(1,1:4)-op_struct.lts_ind(1,1);

length_rx_samples = length(raw_rx_dec)-first_lts_ind-params.TX_NUM_SAMPS+1;%-2*1.5*params.TX_NUM_SAMPS; %max index for the start of a packet
lts_ind_list = first_lts_ind + (0:params.TX_NUM_SAMPS: length_rx_samples);
op_struct.lts_ind_list = lts_ind_list;

% total number of packet and number of packets to process and save
num_packets = length(lts_ind_list);
% params.num_packets = min(num_packets,params.num_packets*params.packet_step_size);
params.num_packets = min(floor(num_packets/params.packet_step_size),params.num_packets);
pack_ind_list = 1:params.packet_step_size:(params.packet_step_size* params.num_packets);

% op_struct.lts_ind_list = lts_ind_list;
% 
% % total number of packet and number of packets to process and save
% num_packets = length(lts_ind_list);
% % params.num_packets = min(num_packets,params.num_packets*params.packet_step_size);
% params.num_packets = num_packets;
% pack_ind_list = 1:params.packet_step_size:(params.packet_step_size* params.num_packets);

op_struct.rx_H_est_vec = zeros(params.num_users,params.NUM_LTS,params.N_SC,params.num_packets);
op_struct.snr = zeros(params.num_users,params.num_packets);
op_struct.packet_payload_list = zeros(params.num_packets,params.num_users);
op_struct.packet_lts_list = zeros(params.num_packets,params.num_users);


op_struct.lts_len = ((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)*params.num_users;
for pidx = 1:num_packets
    %% Correct Freq Offset and Extract Channel

    lts_ind_loop=lts_ind_list(pack_ind_list(pidx));
    try
        % number of original samples ex: 8000 = 80*100
        tx_num_origsamp = (params.N_SC+params.CP_LEN)*params.N_OFDM_SYMS+params.inter_user_zeros;
        %packet index range, ex: rx_lts to rx_lts+160+8000-1
        rx_pkt_and_preamble_range = lts_ind_loop(1):(lts_ind_loop(1)+((params.NUM_LTS+0.5)*(params.N_SC)+params.inter_user_zeros)*params.num_users+tx_num_origsamp-1);
        raw_rx_dec_iter = raw_rx_dec(rx_pkt_and_preamble_range);
%         plot(abs(raw_rx_dec(min(rx_pkt_and_preamble_range)-1000:max(rx_pkt_and_preamble_range)+1000)))
%         xline(1000)
%         xline(max(rx_pkt_and_preamble_range)-(min(rx_pkt_and_preamble_range)-1000))
%         lts_indices_rel_first_user = lts_ind_loop-lts_ind_loop(1)+1;
        [rx_H_est,rx_vec_cfo_corr,params]= extractChannel(raw_rx_dec_iter,1,params,op_struct.user_lts_offsets); 
    catch
        warning('ofdm_rx: Unable to extract channel for a packet');
    end
    op_struct.rx_H_est_vec(:,:,:,pidx) = rx_H_est;
    op_struct.packet_payload_list(pidx)= lts_ind_loop+((params.NUM_LTS+0.5)*(params.N_SC)+params.inter_user_zeros)*params.num_users+1;
    op_struct.packet_lts_list(pidx)= lts_ind_loop;
%     lts_num_users_ind_vec = fliplr(1:1:4);
%     op_struct.packet_payload_list(pidx,:)= lts_ind_loop+((params.NUM_LTS+0.5)*(params.N_SC)+params.inter_user_zeros)*(lts_num_users_ind_vec)+1;
%     op_struct.packet_lts_list(pidx,:)= lts_ind_loop;

%     op_struct.payload_inds =payload_inds;
    
%     op_struct.meansnr(pidx) = params.meansnr; % meansnr obtained by extractchannel()
%     op_struct.CFO_in_Hz(pidx,1) = params.rx_cfo_est_lts*params.SAMP_FREQ;
    
end


%% Plotting
if(params.plot_flag)
    for pidx = 1:10
        timeResol = params.TX_NUM_SAMPS/params.SAMP_FREQ; %sec
        currentTime = timeResol*(pidx-1);
        fprintf('Time index %3.2f ms\n', currentTime*1000)
        figure(params.cf);hold on;
        plot_channel_estimates(op_struct.rx_H_est_vec(:,pidx),params);
        if(params.DO_DECODE)
            figure(params.cf+1);clf;
            plot_evm_and_snr(op_struct.evm_mat(:,:,pidx),op_struct.aevms(pidx),op_struct.snr(pidx),params)
            figure(params.cf+2);clf;
            plot_constellation(params.MOD_ORDER, params.tx_syms_mat,op_struct.rx_syms(:,pidx))
        end
        pause(.5)
    end
end

return

%% Calculate Rx stats
if(params.calc_stats)
    sym_errs = sum(params.tx_data ~= rx_data);
    bit_errs = length(find(dec2bin(bitxor(params.tx_data, rx_data),8) == '1'));
    rx_evm   = sqrt(sum((real(rx_syms) - real(params.tx_syms)).^2 + (imag(rx_syms) - imag(params.tx_syms)).^2)/(length(params.SC_IND_DATA) * params.N_OFDM_SYMS));

    fprintf('\nResults:\n');
    fprintf('Num Bytes:   %d\n', params.N_DATA_SYMS * log2(params.MOD_ORDER) / 8);
    fprintf('Sym Errors:  %d (of %d total symbols)\n', sym_errs, params.N_DATA_SYMS);
    fprintf('Bit Errors:  %d (of %d total bits)\n', bit_errs, params.N_DATA_SYMS * log2(params.MOD_ORDER));

    cfo_est_lts = rx_cfo_est_lts*(params.SAMP_FREQ/params.INTERP_RATE);
    cfo_est_phaseErr = mean(diff(unwrap(pilot_phase_err)))/(4e-6*2*pi);
    cfo_total_ppm = ((cfo_est_lts + cfo_est_phaseErr) /  ((2.412+(.005*(params.CHANNEL-1)))*1e9)) * 1e6;

    fprintf('CFO Est:     %3.2f kHz (%3.2f ppm)\n', (cfo_est_lts + cfo_est_phaseErr)*1e-3, cfo_total_ppm);
    fprintf('     LTS CFO Est:                  %3.2f kHz\n', cfo_est_lts*1e-3);
    fprintf('     Phase Error Residual CFO Est: %3.2f kHz\n', cfo_est_phaseErr*1e-3);

    if params.DO_APPLY_SFO_CORRECTION
        drift_sec = pilot_slope_mat / (2*pi*312500);
        sfo_est_ppm =  1e6*mean((diff(drift_sec) / 4e-6));
        sfo_est = sfo_est_ppm*20;
        fprintf('SFO Est:     %3.2f Hz (%3.2f ppm)\n', sfo_est, sfo_est_ppm);
    end
end

if(params.plot_flag)
    %% Plot Results
%     params.cf= 0;
    % Rx signal
    params.cf= params.cf+ 1;
    figure(params.cf); clf;
    subplot(2,1,1);
    plot(real(rx_vec_air), 'b');
    axis([0 length(rx_vec_air) -params.TX_SCALE params.TX_SCALE])
    grid on;
    title('Rx Waveform (I)');

    subplot(2,1,2);
    plot(imag(rx_vec_air), 'r');
    axis([0 length(rx_vec_air) -params.TX_SCALE params.TX_SCALE])
    grid on;
    title('Rx Waveform (Q)');

    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_rxIQ', params.example_mode_string), '-dpng', '-r96', '-painters')
    end

    % Rx LTS correlation
    params.cf= params.cf+ 1;
    figure(params.cf); clf;
    plot(params.lts_to_plot, '.-b', 'LineWidth', 1);
    hold on;
    grid on;
    line([1 length(params.lts_to_plot)], params.LTS_CORR_THRESH*max(params.lts_to_plot)*[1 1], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2);
    title('LTS Correlation and Threshold')
    xlabel('Sample Index')
    myAxis = axis();
    axis([1, 1000, myAxis(3), myAxis(4)])
    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_ltsCorr', params.example_mode_string), '-dpng', '-r96', '-painters')
    end

    % Channel Estimates
    params.cf= params.cf+ 1;
    figure(params.cf); clf;
    subplot(2,1,1);
    stairs(x - (20/(2*params.N_SC)), fftshift(real(rx_H_est_plot)), 'b', 'LineWidth', 2);
    hold on
    stairs(x - (20/(2*params.N_SC)), fftshift(imag(rx_H_est_plot)), 'r', 'LineWidth', 2);
    hold off
    axis([min(x) max(x) -1.1*max(abs(rx_H_est_plot)) 1.1*max(abs(rx_H_est_plot))])
    grid on;
    title('Channel Estimates (I and Q)')

    subplot(2,1,2);
    bh = bar(x, fftshift(abs(rx_H_est_plot)),1,'LineWidth', 1);
    shading flat
    set(bh,'FaceColor',[0 0 1])
    axis([min(x) max(x) 0 1.1*max(abs(rx_H_est_plot))])
    grid on;
    title('Channel Estimates (Magnitude)')
    xlabel('Baseband Frequency (MHz)')

    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_chanEst', params.example_mode_string), '-dpng', '-r96', '-painters')
    end

    %% Pilot phase error estimate
    params.cf= params.cf+ 1;
    figure(params.cf); clf;
    subplot(2,1,1)
    plot(pilot_phase_err, 'b', 'LineWidth', 2);
    title('Phase Error Estimates')
    xlabel('OFDM Symbol Index')
    ylabel('Radians')
    axis([1 params.N_OFDM_SYMS -3.2 3.2])
    grid on
    h = colorbar;
    set(h,'Visible','off');

    subplot(2,1,2)
    imagesc(1:params.N_OFDM_SYMS, (params.SC_IND_DATA - params.N_SC/2), fftshift(pilot_phase_sfo_corr,1))
    xlabel('OFDM Symbol Index')
    ylabel('Subcarrier Index')
    title('Phase Correction for SFO')
    colorbar
    myAxis = caxis();
    if(myAxis(2)-myAxis(1) < (pi))
       caxis([-pi/2 pi/2])
    end
    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_phaseError', params.example_mode_string), '-dpng', '-r96', '-painters')
    end

    %% Symbol constellation
    params.cf = params.cf+ 1;
    figure(params.cf); clf;
    plot(payload_syms_mat(:),'ro','MarkerSize',1);
    axis square; axis(1.5*[-1 1 -1 1]);
    grid on;
    hold on;
    plot((params.tx_syms_mat(:)+0.00001j),'bo','MarkerSize',5);
    hold off
    title('Tx and Rx Constellations')
    legend('Rx','Tx','Location','EastOutside');

    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_constellations', params.example_mode_string), '-dpng', '-r96', '-painters')
    end

    % EVM & SNR
    params.cf= params.cf+ 1;
    figure(params.cf); clf;
    subplot(2,1,1)
    plot(100*evm_mat(:),'o','MarkerSize',1)
    axis tight
    hold on
    plot([1 length(evm_mat(:))], 100*[aevms, aevms],'r','LineWidth',4)
    myAxis = axis;
    h = text(round(.05*length(evm_mat(:))), 100*aevms+ .1*(myAxis(4)-myAxis(3)), sprintf('Effective SNR: %.1f dB', snr));
    set(h,'Color',[1 0 0])
    set(h,'FontWeight','bold')
    set(h,'FontSize',10)
    set(h,'EdgeColor',[1 0 0])
    set(h,'BackgroundColor',[1 1 1])
    hold off
    xlabel('Data Symbol Index')
    ylabel('EVM (%)');
    legend('Per-Symbol EVM','Average EVM','Location','NorthWest');
    title('EVM vs. Data Symbol Index')
    grid on

    subplot(2,1,2)
    imagesc(1:params.N_OFDM_SYMS, (params.SC_IND_DATA - params.N_SC/2), 100*fftshift(evm_mat,1))
    grid on
    xlabel('OFDM Symbol Index')
    ylabel('Subcarrier Index')
    title('EVM vs. (Subcarrier & OFDM Symbol)')
    h = colorbar;
    set(get(h,'title'),'string','EVM (%)');
    myAxis = caxis();
    if (myAxis(2)-myAxis(1)) < 5
        caxis([myAxis(1), myAxis(1)+5])
    end

    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_evm', params.example_mode_string), '-dpng', '-r96', '-painters')
    end

end

end
