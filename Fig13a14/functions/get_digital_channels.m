function [digital_channel_mat, user_snrs_pre_digital,ret_val] = get_digital_channels(meta_params,ofdm_tx_structs,summed_sigs)
    
    fig_idx = meta_params.fig_idx;
    oversamp_fac = meta_params.oversamp_fac;
    disp_plots = meta_params.plot_setting;
    
    num_users = meta_params.num_users;
    num_ants = meta_params.num_ants;
    num_lts_to_use = 10;
%     disp_plots = true;
    conj_ant = 1;
    snr_calc_ant = 1;
    ofdm_params_to_use = ofdm_tx_structs{1};
    
    clear op_struct_cell op_struct_conj
    
%     for usr_idx=1:1:num_users
    downsamp_conj = summed_sigs{1};    
    [op_struct_conj, ~] = ofdm_rx(downsamp_conj, ofdm_params_to_use);
    if(isempty(op_struct_conj))
        disp("LTS not detected, exiting")
        digital_channel_mat=0;
        user_snrs_pre_digital=0;
        ret_val = -1;
        return;
    end
%     end
    ret_val=1;

    op_struct_samp = op_struct_conj;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct_samp.rx_H_est_vec);
    conj_avg_chans = zeros(num_users,num_users,num_subc,num_packets);
    user_snrs_pre_digital = zeros([num_users,num_users]);
    

    for user_idx=1:1:num_users
        
        curr_conj_struct = op_struct_conj;
        
        lts_inds = curr_conj_struct.lts_ind;
        ref_chan = curr_conj_struct.rx_H_est_vec(:,:,:,start_offset:end);
        ref_mult = exp(-1j*angle(ref_chan));

        [op_struct{user_idx}, ~] = ofdm_rx(summed_sigs{user_idx}, ofdm_params_to_use, lts_inds);

        curr_chan = op_struct{user_idx}.rx_H_est_vec(:,:,:,start_offset:end);
        conj_chan = squeeze(curr_chan.*ref_mult);
%                 conj_chan_lts_avg = squeeze(mean(conj_chan,1));
        conj_chan_lts_avg = squeeze(mean(conj_chan(:,ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:),2));
        conj_avg_chans(user_idx,:,:,:) = conj_chan_lts_avg;


        conj_chan_lts_std = squeeze(std(abs(conj_chan(:,ofdm_params_to_use.NUM_LTS-num_lts_to_use+1:ofdm_params_to_use.NUM_LTS,:,:)),0,2));
        snr_lin_per_subc = abs(conj_chan_lts_avg)./conj_chan_lts_std;
        snr_db_allpackets= 20*log10(mean(snr_lin_per_subc(:,ofdm_params_to_use.SC_IND_DATA,:),[2,3])); % avg across subc and packets
        user_snrs_pre_digital(user_idx,:) = snr_db_allpackets;

%                 chan_pow_allsubc_db = 20*log10(abs(conj_chan_lts_avg));
%                 chan_pow = mean(chan_pow_allsubc_db(ofdm_params_to_use.SC_IND_DATA,:),[1,2]);

        if(disp_plots)
            for user_idx_chan=1:1:num_users
                subplot(num_users,2*num_users,(user_idx-1)*2*num_users+(user_idx_chan-1)*2+1)
                plot(20*log10(abs(squeeze(conj_chan_lts_avg(user_idx_chan,:,:)))),'Color','b')
                title("Relative channel magnitude: RFC: "+num2str(user_idx)+" User idx: "+num2str(user_idx_chan))
                hold on
    %                 ylim([-40,20])
                subplot(num_users,2*num_users,(user_idx-1)*2*num_users+(user_idx_chan)*2)
                plot(unwrap(fftshift(angle(squeeze(conj_chan_lts_avg(user_idx_chan,:,:)))))*180/pi,'Color','b')
                title("Relative channel phase: RFC: "+num2str(user_idx)+" User idx: "+num2str(user_idx_chan))
                ylim([-200,200])
                hold on
            end
        end
    end
    if(disp_plots)
        sgtitle("Channels for user: "+num2str(user_idx))
    end


    
    digital_channel_mat = conj_avg_chans;
end

