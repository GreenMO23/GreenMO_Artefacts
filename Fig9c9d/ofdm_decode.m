function [op_struct, params] = ofdm_decode(rx_vec_air, tx_params, params, user_index, lts_offsets)
%OFDM_RX OFDM SISO receiver, from WARPLab, Generalized to to any number of subcarriers
%   Contributors: Ish Jain, Raghav Subbaraman
    sfo_en = true;
    % total number of packet and number of packets to process and save
    num_packets = 1;
    pidx=1;
    op_struct.rx_H_est_vec = zeros(params.N_SC,num_packets);
    op_struct.rx_syms = zeros(params.N_OFDM_SYMS*params.N_DATA_IND,num_packets);
    op_struct.rx_data = zeros(params.N_OFDM_SYMS*params.N_DATA_IND,num_packets);
    op_struct.evm_mat = zeros([size(params.tx_syms_mat), num_packets]);
    op_struct.aevms = zeros(1,num_packets);
    op_struct.snr = zeros(1,num_packets);
    op_struct.bernum = zeros(1,num_packets);
    op_struct.berratio = zeros(1,num_packets);


    %% Correct Freq Offset and Extract Channel
    % number of original samples ex: 8000 = 80*100
    tx_num_origsamp = (params.N_SC+params.CP_LEN)*params.N_OFDM_SYMS;
    %packet index range, ex: rx_lts to rx_lts+160+8000-1
    rx_pkt_and_preamble_range = 1:(1+((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)*params.num_users+tx_num_origsamp+0.5*params.inter_user_zeros-1);
    raw_rx_dec_iter = rx_vec_air(rx_pkt_and_preamble_range);
    [rx_H_est,rx_vec_cfo_corr,params]= extractChannel_singuser(raw_rx_dec_iter,params,user_index,lts_offsets); 

    op_struct.rx_H_est_vec(:,pidx) = rx_H_est;
    
    op_struct.CFO_in_Hz(pidx,1) = params.rx_cfo_est_lts*params.SAMP_FREQ;

    %% Rx Payload
    syms_eq_mat = rx_equalize_signal(rx_vec_cfo_corr,rx_H_est,params,lts_offsets,user_index);
    % SFO and phase correction
    if(sfo_en)
        [syms_eq_pc_mat,params] = sfo_and_phase_error_correction(syms_eq_mat,params);

         %% Demodulate 
         payload_syms_mat = syms_eq_pc_mat(params.SC_IND_DATA, :);
    else
        payload_syms_mat = syms_eq_mat(params.SC_IND_DATA, :);
    end

    if(~params.enable_channel_codes)
        op_struct.rx_syms(:,pidx) = reshape(payload_syms_mat, params.N_DATA_SYMS,1);
        % qamdemod is a very slow function. Consider replacing with faster
        % alternative - Raghav
        op_struct.rx_data(:,pidx) = qamdemod(op_struct.rx_syms(:,pidx),params.MOD_ORDER,'UnitAveragePower',true);
        tx_syms_mat = tx_params.tx_syms_mat;
        tx_data = tx_params.tx_data;
        %% Calculate Stats
        op_struct.evm_mat(:,:,pidx) = abs(payload_syms_mat - tx_syms_mat).^2;
        op_struct.aevms(pidx)= mean(op_struct.evm_mat(:,:,pidx),'all');
    %     op_struct.snr(pidx) = 10*log10(1./op_struct.aevms(pidx)); % dont do
    %     log
        op_struct.snr(pidx) = 1./op_struct.aevms(pidx);

        mod_bit_size = log2(params.MOD_ORDER);
        tx_bits = de2bi(tx_data,mod_bit_size);
        rx_bits = de2bi(op_struct.rx_data(:,pidx),mod_bit_size);
        [op_struct.bernum(pidx),op_struct.berratio(pidx)] = biterr(tx_bits,rx_bits);
    else
        tx_syms_mat = tx_params.tx_syms_mat;
%         tx_data = tx_params.tx_data;
        % Calculate Stats
        op_struct.evm_mat(:,:,pidx) = abs(payload_syms_mat - tx_syms_mat).^2;
        op_struct.aevms(pidx)= mean(op_struct.evm_mat(:,:,pidx),'all');
    %     op_struct.snr(pidx) = 10*log10(1./op_struct.aevms(pidx)); % dont do
    %     logrx
        op_struct.snr(pidx) = 1./op_struct.aevms(pidx);
        
        op_struct.rx_syms(:,pidx) = reshape(payload_syms_mat, params.N_DATA_SYMS,1);
        % qamdemod is a very slow function. Consider replacing with faster
        % alternative - Raghav
        rx_bits_raw = qamdemod(op_struct.rx_syms(:,pidx),params.MOD_ORDER, 'OutputType','bit','UnitAveragePower',true);
%         mod_bit_size = log2(params.MOD_ORDER);
%         tx_bits = de2bi(tx_data,mod_bit_size);
%         rx_bits_raw = de2bi(op_struct.rx_data(:,pidx),mod_bit_size).';
        rx_bits_flattened = rx_bits_raw(:).';
%         biterr(tx_params.tx_data,rx_bits_flattened);
        
        if(tx_params.channel_coding_rate==0.5)
            % viterbi decoder
            trel = poly2trellis(7, [171 133]);              % Define trellis
            rx_bits = vitdec(rx_bits_flattened,trel,7,'trunc','hard');
            tx_bits = tx_params.tx_bits;
            [err_num, err_rat] = biterr(tx_bits,rx_bits(1:numel(tx_bits)));
            op_struct.bernum(pidx) = err_num;
            op_struct.berratio(pidx) = err_rat;
        else
            vitDecoder = comm.ViterbiDecoder(poly2trellis(7, [171 133]), ...
            'InputFormat', 'Hard');
            vitDecoder.PuncturePatternSource =  'Property';
            vitDecoder.PuncturePattern = [1;1;0;1;1;0];
            vitDecoder.TracebackDepth = 96;
            errorCalc = comm.ErrorRate('ReceiveDelay', vitDecoder.TracebackDepth);
            decData = vitDecoder(rx_bits_flattened.');
            % Compute and accumulate errors
            tx_bits = tx_params.tx_bits;
            err_stat = errorCalc(tx_bits.', decData);
            err_rat = err_stat(1);
            err_num = err_stat(2);
            op_struct.bernum(pidx) = err_num;
            op_struct.berratio(pidx) = err_rat;
            reset(vitDecoder);
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
    grid on;;
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
