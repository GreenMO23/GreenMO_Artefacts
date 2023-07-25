function [code_shifted_rx_sigs] = apply_delayed_clocks_twouser(meta_params,rx_interfered_waveforms,channel_est_struct)
    oversamp_fac = meta_params.oversamp_fac;
    [code_len,~] = size(rx_interfered_waveforms{1});
    
    shift_freq = 1; %1x sampling freq
    duty = 0.5;
    code_num_samples = floor(oversamp_fac/shift_freq);
    code_vec_1r = [ones([1,duty*code_num_samples]) zeros([1,(1-duty)*code_num_samples])];
    num_code_reps = floor(code_len/code_num_samples);
    code_vec = repmat(code_vec_1r,[1,num_code_reps]);
    base_clk_rep{1} = code_vec;

    shift_freq = 2; %1x sampling freq
    duty = 0.5;
    code_num_samples = floor(oversamp_fac/shift_freq);
    code_vec_1r_2 = [ones([1,duty*code_num_samples]) zeros([1,(1-duty)*code_num_samples])];
    num_code_reps = floor(code_len/code_num_samples);
    code_vec_2 = repmat(code_vec_1r_2,[1,num_code_reps]);
    base_clk_rep{2} = code_vec_2;
    
    num_ants = channel_est_struct.num_ants;
    see_spectrum = false;
    ant_idx_to_see = 3;
    if(mod(code_len,oversamp_fac)==0)
        
        clear code_shifted_rx_sigs code_pattern rx_interfered_waveforms_downsampled
        for ant_idx = 1:1:num_ants
            for usr_idx=1:1:meta_params.num_users
                curr_delays = channel_est_struct.delays{usr_idx};
                if(usr_idx==1)
                    code_pattern{ant_idx} = delayseq2(base_clk_rep{usr_idx}.',curr_delays(ant_idx));
                else
                    code_pattern{ant_idx} = code_pattern{ant_idx}+ delayseq2(base_clk_rep{usr_idx}.',curr_delays(ant_idx));
                end
                
            end
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
                sgtitle("Codes for Antenna: "+num2str(ant_idx))
                
                figure(11)
                subplot(2,1,1)
                abs_vec = abs(fftshift(fft(code_pattern{ant_idx})));
                mask_vec = abs_vec>0.1*max(abs_vec);
                plot(freq_vec,abs_vec)
                subplot(2,1,2)
                angle_vec = angle(fftshift(fft(code_pattern{ant_idx}))).*mask_vec;
                plot(freq_vec,angle_vec*180/pi)
                xlim([-2.5,2.5])
                ylim([-180,180])
                print("debug");
                
%                 subplot(2,2,3)
%                 abs_vec = abs(fftshift(fft(base_clk_rep{2})));
%                 mask_vec = abs_vec>0.1*max(abs_vec);
%                 plot(freq_vec,abs_vec)
%                 subplot(2,2,4)
%                 angle_vec = angle(fftshift(fft(base_clk_rep{2}))).*mask_vec;
%                 plot(freq_vec,angle_vec)
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

