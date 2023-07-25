function [code_shifted_rx_sigs] = apply_delayed_clocks(meta_params,rx_interfered_waveforms,channel_est_struct)
    oversamp_fac = meta_params.oversamp_fac;
    [code_len,~] = size(rx_interfered_waveforms{1});
    num_ants = channel_est_struct.num_ants;
    base_clk = zeros([oversamp_fac,1]);
    base_clk(1) = 1;
    see_spectrum = false;
    ant_idx_to_see = 1;
    if(mod(code_len,oversamp_fac)==0)
        base_clk_rep = repmat(base_clk,[floor(code_len/oversamp_fac),1]);
        clear code_shifted_rx_sigs code_pattern
        for ant_idx = 1:1:num_ants
            code_pattern{ant_idx} = delayseq2(base_clk_rep,channel_est_struct.delays(ant_idx));
            code_shifted_rx_sigs{ant_idx} = rx_interfered_waveforms{ant_idx}.*code_pattern{ant_idx};
            if(see_spectrum && ant_idx==ant_idx_to_see)
                figure(10)
                freq_vec = linspace(-oversamp_fac/2,oversamp_fac/2,code_len);
                subplot(2,2,1)
                plot(freq_vec,abs(fftshift(fft(rx_interfered_waveforms{ant_idx}))))
                subplot(2,2,2)
                plot(freq_vec,angle(fftshift(fft(rx_interfered_waveforms{ant_idx}))))
                subplot(2,2,3)
                plot(freq_vec,abs(fftshift(fft(code_shifted_rx_sigs{ant_idx}))))
                subplot(2,2,4)
                plot(freq_vec,angle(fftshift(fft(code_shifted_rx_sigs{ant_idx}))))
                sgtitle("Channels for Antenna: "+num2str(ant_idx))
                
                figure(11)
                subplot(2,2,1)
                abs_vec = abs(fftshift(fft(base_clk_rep)));
                mask_vec = abs_vec>0.1*max(abs_vec);
                plot(freq_vec,abs_vec)
                subplot(2,2,2)
                angle_vec = angle(fftshift(fft(base_clk_rep))).*mask_vec;
                plot(freq_vec,angle_vec)
                
                delayed_codes = delayseq(base_clk_rep,pi/2*(oversamp_fac/2*pi));
                subplot(2,2,3)
                abs_vec = abs(fftshift(fft(delayed_codes)));
                mask_vec = abs_vec>0.1*max(abs_vec);
                plot(freq_vec,abs_vec)
                subplot(2,2,4)
                angle_vec = angle(fftshift(fft(delayed_codes))).*mask_vec;
                plot(freq_vec,angle_vec)
            end
        end
    else
        disp("Not implemented for this usecase")
        return
    end
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins.'));
end

