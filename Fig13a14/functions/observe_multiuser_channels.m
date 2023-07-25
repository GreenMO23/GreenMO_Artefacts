function [digital_channel_mat, user_snrs_pre_digital] = observe_multiuser_channels(meta_params,ofdm_tx_structs,channel_struct,resampled_windowed_user_sigs);
    
    fig_idx = meta_params.fig_idx;
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
    
%     for usr_idx=1:1:num_users
    downsamp_conj = resampled_windowed_user_sigs{1,conj_ant};    
    [op_struct_conj, ~] = ofdm_rx(downsamp_conj, ofdm_params_to_use);
    if(isempty(op_struct_conj))
        disp("LTS not detected, exiting")
        chan_est = -1;
        return;
    end
%     end
    

    op_struct_samp = op_struct_conj;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct_samp.rx_H_est_vec);
    conj_avg_chans = zeros(num_ants,num_users,num_subc,num_packets);
    user_snrs_pre_digital = zeros([num_users,num_users]);
    

    for user_idx=1:1:num_users
        if(disp_plots)
            figure(user_idx+fig_idx)
        end
        curr_conj_struct = op_struct_conj;
        for ant_idx=1:1:num_ants
            lts_inds = curr_conj_struct.lts_ind;
            ref_chan = curr_conj_struct.rx_H_est_vec(:,:,:,start_offset:end);
            ref_mult = exp(-1j*angle(ref_chan));

            [op_struct{user_idx,ant_idx}, ~] = ofdm_rx(resampled_windowed_user_sigs{user_idx,ant_idx}, ofdm_params_to_use, lts_inds);
            
            curr_chan = op_struct{user_idx,ant_idx}.rx_H_est_vec(:,:,:,start_offset:end);
            conj_chan = squeeze(curr_chan.*ref_mult);
%                 conj_chan_lts_avg = squeeze(mean(conj_chan,1));
            conj_chan_lts_avg = squeeze(mean(conj_chan(:,ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:),2));
            conj_avg_chans(ant_idx,:,:,:) = conj_chan_lts_avg;

            if(ant_idx==snr_calc_ant)
                
                conj_chan_lts_std = squeeze(std(abs(conj_chan(:,ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:)),0,2));
                snr_lin_per_subc = abs(conj_chan_lts_avg)./conj_chan_lts_std;
                snr_db_allpackets= 20*log10(mean(snr_lin_per_subc(:,ofdm_params_to_use.SC_IND_DATA,:),[2,3])); % avg across subc and packets
                user_snrs_pre_digital(user_idx,:) = snr_db_allpackets;
                
%                 chan_pow_allsubc_db = 20*log10(abs(conj_chan_lts_avg));
%                 chan_pow = mean(chan_pow_allsubc_db(ofdm_params_to_use.SC_IND_DATA,:),[1,2]);
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
    
    digital_channel_mat = conj_avg_chans;
    
end

