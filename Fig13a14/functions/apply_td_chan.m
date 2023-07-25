function [rx_waveforms_interfered, rx_waveforms_interfered_no_noise,rx_waveforms_interfered_dig, rx_waveforms_interfered_babf] = apply_td_chan(meta_params,tx_waveforms, loopback,channel_struct)
    clear rx_waveforms_interfered rx_waveforms_interfered_dig rx_waveforms_interfered_babf rx_waveforms_interfered_no_noise
    
%     noise_pow = 0.001*10^(meta_params.tx_pow-meta_params.user1_SNR);
    
    for i = 1:1:channel_struct.num_users
%         disp("Applying channel for user: "+num2str(i));
        for j=1:1:channel_struct.num_channel_taps
            for k=1:1:channel_struct.num_ants
                delay_val = channel_struct.taps_delays(i,j,k)*meta_params.bandwidth*meta_params.oversamp_fac*(10^-9);% taps_delays is in ns
                
%                 if(delay_val>channel_struct.max_delay)
%                     disp("Not designed for such high delays, exiting")
%                     return 
%                 end
                oversamp_tx = resample(tx_waveforms{i},meta_params.oversamp_fac,1)*channel_struct.taps_phases(i,j,k);
%                 oversamp_loopback = resample(loopback{i},meta_params.oversamp_fac,1)*channel_struct.taps_phases(i,j,k);
%                 oversamp_tx = interp(tx_waveforms{i},meta_params.oversamp_fac)*channel_struct.taps_phases(i,j,k);
%                 oversamp_tx = tx_waveforms{i}*channel_struct.taps_phases(i,j,k);
                [num_samps,~] = size(oversamp_tx);
                td_sig = delayseq2(oversamp_tx*channel_struct.taps_mag(i,j,k),delay_val);
%                 if(j==1)
%                     rx_waveforms_per_user{i,k}=delayseq2(oversamp_loopback *channel_struct.taps_mag(i,j,k),delay_val);
%                 else
%                     rx_waveforms_per_user{i,k}=rx_waveforms_per_user{i,k}+delayseq2(oversamp_loopback*channel_struct.taps_mag(i,j,k),delay_val);
%                 end
                
                if(i==1 && j==1) % first tap and first user
%                     noise_sig = sqrt(noise_pow/2)*(randn(num_samps,1)+1j*randn(num_samps,1));
                    rx_waveforms_interfered{k} = td_sig;
                    rx_waveforms_interfered_no_noise{k} = td_sig;
                    rx_waveforms_interfered_dig{k} = td_sig;
                    rx_waveforms_interfered_babf{k} = td_sig;
                else
                    % add for every other tap and every other user
                    rx_waveforms_interfered{k} = rx_waveforms_interfered{k}+td_sig; 
                    rx_waveforms_interfered_no_noise{k} = rx_waveforms_interfered_no_noise{k}+td_sig; 
                    rx_waveforms_interfered_dig{k} = rx_waveforms_interfered_dig{k}+td_sig;
                    rx_waveforms_interfered_babf{k} = rx_waveforms_interfered_babf{k}+td_sig;
%                     rx_waveforms_interfered_dig{k} = rx_waveforms_interfered_dig{k}+delayseq2(oversamp_tx*channel_struct.taps_mag(i,j,k),delay_val); 
%                     rx_waveforms_interfered_babf{k} = rx_waveforms_interfered_babf{k}+delayseq2(oversamp_tx*channel_struct.taps_mag(i,j,k),delay_val); 
                end
                
                % add noise accoridngly
                if(i==1 && j==channel_struct.num_channel_taps)
                    
                    first_ltses = rx_waveforms_interfered{k};
                    lts_sensitivity = 10;
                    first_ltses = first_ltses(501+lts_sensitivity:500+10.5*64*meta_params.oversamp_fac-lts_sensitivity);
                    rx_pow = rms(first_ltses)^2;
                    noise_pow = rx_pow*(10^((-meta_params.user1_SNR)/10));
                    noise_sig = sqrt(noise_pow/2)*(randn(num_samps,1)+1j*randn(num_samps,1));
%                     base_sig = rx_waveforms_interfered{k};
                    rx_waveforms_interfered_dig{k} = rx_waveforms_interfered_dig{k}+noise_sig*sqrt(channel_struct.num_users);
                    rx_waveforms_interfered_babf{k} = rx_waveforms_interfered_babf{k}+noise_sig;
                    rx_waveforms_interfered{k} = rx_waveforms_interfered{k}+noise_sig;
%                     rx_waveforms_per_user{i,k} = base_sig+noise_sig; 
%                     figure(44)
%                     plot(abs(rx_waveforms_interfered{k}))
%                     hold on
%                     xline(501+lts_sensitivity);
%                     xline(500+10.5*64*meta_params.oversamp_fac-lts_sensitivity);
%                     title("Debug noise add");

%                     rx_waveforms_interfered{k} = sinc_interp(rx_waveforms_interfered{k},meta_params.oversamp_fac);                    
%                     first_ltses = rx_waveforms_interfered{k};
%                     lts_sensitivity = 10;
%                     first_ltses = first_ltses(1+lts_sensitivity:10.5*64*meta_params.oversamp_fac-lts_sensitivity);
%                     curr_sig = first_ltses;
%                     curr_sample_size = numel(curr_sig);
%                     num_zeros_remove = mod(curr_sample_size,64*meta_params.oversamp_fac);
%                     curr_sig_zero_removed = curr_sig(1:curr_sample_size-num_zeros_remove);
%                     curr_sig_reshaped = reshape(curr_sig_zero_removed,64*meta_params.oversamp_fac,[]);
%                     fft_mat = fft(curr_sig_reshaped,64*meta_params.oversamp_fac,1);
%                     rx_pow = mean(abs(fft_mat(1:32,1));
%                     noise_pow = 10^((rx_pow-meta_params.user1_SNR)/10);
%                     noise_sig = sqrt(noise_pow/2)*(randn(num_samps,1)+1j*randn(num_samps,1));
%                     plot(20*log10(abs(fft_mat)))
                end
                
            end
        end
    end
end


function interp_sig = sinc_interp(curr_sig,oversamp_fac)
    nfft = numel(curr_sig);
    oversamped_nfft = oversamp_fac*nfft;
    curr_fft = fft(curr_sig);
    oversamped_fft=zeros([oversamped_nfft,1]);
    oversamped_fft(1:nfft/2) = curr_fft(1:nfft/2);
    oversamped_fft(oversamped_nfft-nfft/2+1:oversamped_nfft) = curr_fft(nfft/2+1:nfft);
    interp_sig = ifft(oversamped_fft)*oversamp_fac;
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    if(frac_del>1)
        num_zeros_prepend = floor(frac_del);
        frac_del = mod(frac_del,1);
    else
        num_zeros_prepend = 0;
    end
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig_nozeros = ifft(X .* exp(-1j * frac_del * fbins));
    if(num_zeros_prepend==0)
       delayed_sig = delayed_sig_nozeros.';
    else
        delayed_sig = [zeros(1,num_zeros_prepend) delayed_sig_nozeros(1:nfft-num_zeros_prepend)].';
    end
end


