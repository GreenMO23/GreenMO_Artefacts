function [snr_db,chan_pow] = get_increased_snr(meta_params,ofdm_tx_structs,channel_struct,summed_sigs)
    
    downsample_flag = meta_params.downsamp_flag;
    oversamp_fac = meta_params.oversamp_fac;
    fig_idx = meta_params.fig_idx;
    
    
    num_users = channel_struct.num_users;
    num_ants = channel_struct.num_ants;
    num_lts_to_use = 10;
    disp_plots = meta_params.plot_setting;
    conj_ant = 1;
    ofdm_params_to_use = ofdm_tx_structs{1};
    
    clear op_struct_cell op_struct_conj
    if(downsample_flag==0)
        downsamp_conj = resample(summed_sigs,1,oversamp_fac);
    else
        downsamp_conj = downsamp_filt(summed_sigs,oversamp_fac);
    end
    [op_struct_conj, ~] = ofdm_rx(downsamp_conj, ofdm_params_to_use);

    if(isempty(op_struct_conj))
        disp("LTS not detected, exiting")
        chan_est = -1;
        return;
    end
    
    lts_inds = op_struct_conj.lts_ind;
    op_struct = op_struct_conj;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct.rx_H_est_vec);
    conj_avg_chans = zeros(num_packets,num_users,num_ants,num_subc);
%     conj_idx =3;
    

    curr_chan = op_struct_conj.rx_H_est_vec(1,:,:,start_offset:end);
    conj_chan = squeeze(curr_chan);
%                 conj_chan_lts_avg = squeeze(mean(conj_chan,1));
    conj_chan_lts_avg = squeeze(mean(conj_chan(ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:),1));

    
%     curr_sigs = downsamp_conj;
%     extraction_sensitivity = 10;
%     extracted_lts = curr_sigs(op_struct_conj.lts_ind(1)-64*10+extraction_sensitivity:op_struct_conj.lts_ind(1)+32-extraction_sensitivity);
%     extracted_noise = curr_sigs(op_struct_conj.lts_ind(1)+32+extraction_sensitivity:op_struct_conj.lts_ind(1)+32+ofdm_params_to_use.inter_user_zeros-extraction_sensitivity);
% %     snr_db = 20*log10((rms(extracted_lts)-rms(extracted_noise))/rms(extracted_noise));
%     snr_db = 10*log10((rms(extracted_lts)^2-rms(extracted_noise)^2)/rms(extracted_noise)^2);
    
    conj_chan_lts_std = squeeze(std(abs(conj_chan(ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:)),1));
    snr_lin_per_subc = abs(conj_chan_lts_avg)./conj_chan_lts_std;
    snr_db_allpackets= 20*log10(mean(snr_lin_per_subc(ofdm_params_to_use.SC_IND_DATA,:)));
    snr_db = mean(snr_db_allpackets);

    chan_pow_allsubc_db = 20*log10(abs(conj_chan_lts_avg));
    chan_pow = mean(chan_pow_allsubc_db(ofdm_params_to_use.SC_IND_DATA,:),[1,2]);
    
    
    if(disp_plots)
        figure(fig_idx)
        subplot(1,2,1)
        plot(20*log10(abs(conj_chan_lts_avg)),'Color','b')
        title("Relative channel magnitude: Antenna")
        hold on
%                 ylim([-40,20])
        subplot(1,2,2)
        plot(unwrap(fftshift(angle(conj_chan_lts_avg(ofdm_params_to_use.SC_IND_DATA,:))))*180/pi,'Color','b')
        title("Relative channel phase: Antenna")
%                 ylim([-200,200])
        hold on
    end
    
    


end

function [downsamples] = downsamp_filt(samples,oversamp_fac)
    plot_debug=false;
    orig_fft = fft(samples);
    [num_over_samps,~] = size(samples);
    num_downsamps = num_over_samps/oversamp_fac;
    selected_freq_samples = num_downsamps/(oversamp_fac):num_downsamps/(oversamp_fac)+num_downsamps-1;
    downsamples = ifft(ifftshift(orig_fft(selected_freq_samples)));
    if(plot_debug)
        over_freq_vec = linspace(-oversamp_fac*0.5,oversamp_fac*0.5,num_over_samps);
        subplot(2,1,1)
        plot(abs(orig_fft))
        hold on
        xline(num_downsamps/(oversamp_fac))
        xline(num_downsamps/(oversamp_fac)+num_downsamps-1)
        subplot(2,1,2)
        plot(angle(orig_fft))
    end
end

