function [rx_interfered_waveforms_downsampled] = get_downsampled_wavs(meta_params,rx_interfered_waveforms)
    num_users = meta_params.num_users;
    num_ants = meta_params.num_ants;
    oversamp_fac = meta_params.oversamp_fac;
    
    clear rx_interfered_waveforms_downsampled;
    
    for ant_idx=1:1:num_ants
        rx_interfered_waveforms_downsampled{ant_idx} = resample(rx_interfered_waveforms{ant_idx},1,oversamp_fac);
    end

end

