function [fully_connected_rx_sigs,partially_connected_rx_sigs] = hybrid_beamf_combiner(meta_params,rx_interfered_waveforms_downsampled,channel_est_struct)  
    num_ants = channel_est_struct.num_ants;
    clear fully_connected_rx_sigs partially_connected_rx_sigs
    
    partial_ant_vec=1:1:channel_est_struct.num_ants;
    hyb_num_rfc = 8;
    num_ants_per_rfc = channel_est_struct.num_ants/8;
    partial_ant_vec = reshape(partial_ant_vec,[num_ants_per_rfc,hyb_num_rfc]).';
    
    for usr_idx=1:1:meta_params.num_users
        fully_connected_rx_sigs{usr_idx} = 0;
        curr_beamf_vec = squeeze(channel_est_struct.fully_connected_mat(:,usr_idx));
        for ant_idx = 1:1:num_ants    
            phase_mult = curr_beamf_vec(ant_idx);
            fully_connected_rx_sigs{usr_idx} = fully_connected_rx_sigs{usr_idx}+rx_interfered_waveforms_downsampled{ant_idx}*phase_mult;
        end
        
        partially_connected_rx_sigs{usr_idx}=0;
        for ant_idx = partial_ant_vec(usr_idx,:)
            hyb_phase_mult = curr_beamf_vec(ant_idx);
            partially_connected_rx_sigs{usr_idx} = partially_connected_rx_sigs{usr_idx}+rx_interfered_waveforms_downsampled{ant_idx}*hyb_phase_mult;
        end
%         partially_connected_rx_sigs=-1;
    end
end

