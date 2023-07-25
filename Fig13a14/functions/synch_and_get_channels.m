function [chan_est_struct,chan_est, snr_db] = synch_and_get_channels(meta_params,ofdm_tx_structs,channel_struct,rx_interfered_waveforms)
    
    fig_idx = meta_params.fig_idx;
    downsample_flag = meta_params.downsamp_flag;
    oversamp_fac = meta_params.oversamp_fac;
    
    disp_plots = meta_params.plot_setting;
    
    num_users = channel_struct.num_users;
    num_ants = channel_struct.num_ants;
    num_lts_to_use = 10;
%     disp_plots = true;
    conj_ant = 1;
    snr_calc_ant = 1;
    ofdm_params_to_use = ofdm_tx_structs{1};
    
    clear op_struct_cell op_struct_conj
    
    downsamp_conj = resample(rx_interfered_waveforms{conj_ant},1,oversamp_fac);
    
    [op_struct_conj, ~] = ofdm_rx(downsamp_conj, ofdm_params_to_use);

    if(isempty(op_struct_conj))
        disp("LTS not detected, exiting")
        chan_est = -1;
        return;
    end
    
    lts_inds = op_struct_conj.lts_ind;
    
    clear downsampled_sigs 
    for ant_ind=1:1:num_ants
        downsamp_samp = resample(rx_interfered_waveforms{ant_ind},1,oversamp_fac);
        downsampled_sigs{ant_ind} = downsamp_samp;
        [op_struct_cell{ant_ind}, ~] = ofdm_rx(downsamp_samp, ofdm_params_to_use, lts_inds);
    end

    op_struct = op_struct_cell;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct{1}.rx_H_est_vec);
    conj_avg_chans = zeros(num_packets,num_users,num_ants,num_subc);
%     conj_idx =3;
    
    ref_chan = op_struct{conj_ant}.rx_H_est_vec(1,:,:,start_offset:end);
    ref_mult = exp(-1j*angle(ref_chan));
        
    for user_idx=1:1:num_users
        if(disp_plots)
            figure(user_idx+3+fig_idx)
        end
        

        for ant_idx=1:1:num_ants
            curr_chan = op_struct{ant_idx}.rx_H_est_vec(user_idx,:,:,start_offset:end);
            conj_chan = squeeze(curr_chan.*ref_mult);
%                 conj_chan_lts_avg = squeeze(mean(conj_chan,1));
            conj_chan_lts_avg = squeeze(mean(conj_chan(ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:),1));
            conj_avg_chans(:,user_idx,ant_idx,:) = conj_chan_lts_avg.';
%             if(ant_idx==1)
%                 conj_chan_lts_std = squeeze(std(conj_chan(ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:),1));
%                 snr_lin_per_subc = abs(conj_chan_lts_avg)./conj_chan_lts_std;
%                 snr_db_allpackets= 20*log10(mean(snr_lin_per_subc(ofdm_params_to_use.SC_IND_DATA,:)));
%                 snr_db = mean(snr_db_allpackets);
%             end
            curr_sigs = downsampled_sigs{ant_idx};
%             figure(99)
%             plot(abs(curr_sigs))
%             hold on
%             xline(op_struct{ant_idx}.lts_ind(1)-64*10)
%             xline(op_struct{ant_idx}.lts_ind(1)+32)
%             xline(op_struct{ant_idx}.lts_ind(1)+32+ofdm_params_to_use.inter_user_zeros)
            if(ant_idx==snr_calc_ant)
%                 extraction_sensitivity = 10;
%                 extracted_lts = curr_sigs(op_struct{ant_idx}.lts_ind(1)-64*10+extraction_sensitivity:op_struct{ant_idx}.lts_ind(1)+32-extraction_sensitivity);
%                 extracted_noise = curr_sigs(op_struct{ant_idx}.lts_ind(1)+32+extraction_sensitivity:op_struct{ant_idx}.lts_ind(1)+32+ofdm_params_to_use.inter_user_zeros-extraction_sensitivity);
%                 snr_db = 10*log10((rms(extracted_lts)^2-rms(extracted_noise)^2)/rms(extracted_noise)^2);
                
                conj_chan_lts_std = squeeze(std(abs(conj_chan(ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:)),0,1));
                snr_lin_per_subc = abs(conj_chan_lts_avg)./conj_chan_lts_std;
                snr_db_allpackets= 20*log10(mean(snr_lin_per_subc(ofdm_params_to_use.SC_IND_DATA,:)));
                snr_db(user_idx) = mean(snr_db_allpackets);
                
                chan_pow_allsubc_db = 20*log10(abs(conj_chan_lts_avg));
                chan_pow = mean(chan_pow_allsubc_db(ofdm_params_to_use.SC_IND_DATA,:),[1,2]);
            end
            
            if(disp_plots)
                subplot(num_ants,2,(ant_idx-1)*2+1)
                plot(20*log10(abs(conj_chan_lts_avg)),'Color','b')
                title("Relative channel magnitude: Antenna: "+num2str(ant_idx-1))
                hold on
%                 ylim([-40,20])
                subplot(num_ants,2,(ant_idx-1)*2+2)
                plot(unwrap(fftshift(angle(conj_chan_lts_avg(ofdm_params_to_use.SC_IND_DATA,:))))*180/pi,'Color','b')
                title("Relative channel phase: Antenna: "+num2str(ant_idx-1))
%                 ylim([-200,200])
                hold on
            end
        end
        if(disp_plots)
            sgtitle("Channels for user: "+num2str(user_idx))
        end
    end
    
    
    chan_est = conj_avg_chans;
    delay_res_samples = meta_params.delay_resolution*1e-9*meta_params.bandwidth;
    chan_est_struct.phase_angles = mean(unwrap(angle(squeeze(conj(chan_est(1,1,:,ofdm_params_to_use.SC_IND_DATA))))),2);
%   chan_est_struct.delays = (chan_est_struct.phase_angles/(2*pi*channel_struct.freq))*channel_struct.samp_freq;
    chan_est_struct.ideal_delays = (chan_est_struct.phase_angles/(2*pi))*oversamp_fac;
    if(meta_params.delay_resolution==0)
        chan_est_struct.delays = chan_est_struct.ideal_delays;
    else
        chan_est_struct.delays = chan_est_struct.ideal_delays-mod(chan_est_struct.ideal_delays,delay_res_samples);
    end
    chan_est_struct.num_lts_to_use = num_lts_to_use;
    chan_est_struct.num_ants = num_ants;


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

