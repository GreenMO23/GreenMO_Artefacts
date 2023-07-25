function [resampled_windowed_user_sigs] = reshape_window_downsample(meta_params,code_phase_shifted_rx_sigs)
    num_users = meta_params.num_users;
    num_ants = meta_params.num_ants;
    oversamp_fac = meta_params.oversamp_fac;
    plot_debug = false;
    clear resampled_windowed_user_sigs;
    for usr_idx=1:1:num_users
        for ant_idx=1:1:num_ants
            samples = code_phase_shifted_rx_sigs{ant_idx};
            orig_fft = fft(samples);
            [num_over_samps,~] = size(samples);
            num_downsamps = num_over_samps/oversamp_fac;
            selected_freq_samples = (num_downsamps)*(usr_idx-1)+0.5*num_downsamps:(num_downsamps)*(usr_idx-1)+1.5*num_downsamps-1;
            if(usr_idx==1)
                downsamples = ifft(fftshift(orig_fft(selected_freq_samples)));
            else
                downsamples = ifft(ifftshift(orig_fft(selected_freq_samples)));
            end
            resampled_windowed_user_sigs{usr_idx,ant_idx} = downsamples;
            if(plot_debug && ant_idx==2)
                if(usr_idx==1)
                    figure(99)
                    over_freq_vec = linspace(-oversamp_fac*0.5,oversamp_fac*0.5,num_over_samps);
                    subplot(2,1,1)
                    plot(abs(orig_fft))
                     hold on
                    subplot(2,1,2)
                    plot(angle(orig_fft))
                end
                subplot(2,1,1)
                xline((num_downsamps)*(usr_idx-1)+0.5*num_downsamps)
                xline((num_downsamps)*(usr_idx-1)+1.5*num_downsamps-1)

            end
        end
    end
end

