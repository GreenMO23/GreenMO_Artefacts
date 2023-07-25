function [kmeans_rate_samples] = kmeans_switcher_multirate(meta_params,rx_interfered_waveforms_babf,channel_est_struct,rate_fac)
    num_rfc = rate_fac;
    num_users = meta_params.num_users;
    
    oversamp_fac = meta_params.oversamp_fac/(num_rfc*num_users);
    
    disp_plots = meta_params.plot_setting;
    
    
    num_ants = meta_params.num_ants;
    
    curr_chan = channel_est_struct.kmeans_chan_mat;
    
    
    
    switch_mat = zeros([num_users,rate_fac,num_ants]);
    phase_vec = zeros([num_users,rate_fac]);
    
%     num_classes_kmeans = num_rfc+1; % ignore one group
%     for user_idx=1:1:num_users
%         curr_usr_chan = squeeze(curr_chan(user_idx,:,:));
%         [kmeans_class,kmeans_centers] = kmeans(angle(curr_usr_chan),num_classes_kmeans);
%         class_tally = accumarray(sort(kmeans_class),1);
%         class_vec = 1:1:num_classes_kmeans;
%         [~,min_tally_class] = min(class_tally);
%         class_vec(min_tally_class)=[];
%         for class_idx = 1:1:num_rfc
%             class_val = class_vec(class_idx);
%             switch_mat(user_idx,class_idx,:) = kmeans_class==class_val;
%             phase_vec(user_idx,class_idx) = kmeans_centers(class_val,1);
%         end  
%     end
    num_classes_kmeans = num_rfc; % ignore one group
    for user_idx=1:1:num_users
        curr_usr_chan = squeeze(curr_chan(user_idx,:,:));
        [kmeans_class,kmeans_centers] = kmeans(angle(curr_usr_chan),num_classes_kmeans);
        for class_val = 1:1:num_rfc
            switch_mat(user_idx,class_val,:) = kmeans_class==class_val;
%             phase_vec(user_idx,class_val) = kmeans_centers(class_val,1);
            if(class_val==1)
                curr_chan_summed = sum(squeeze(curr_usr_chan(:,1)).*(kmeans_class==class_val));
                class1_ang = angle(curr_chan_summed);
                phase_vec(user_idx,class_val) = 0;
            else
                curr_chan_summed = sum(squeeze(curr_usr_chan(:,1)).*(kmeans_class==class_val));
                phase_vec(user_idx,class_val) = class1_ang-angle(curr_chan_summed);
            end
        end  
    end
    
%     
    for usr_idx=1:1:num_users
        phase_combined_sigs = 0;
        for rfc_idx=1:1:num_rfc
            summed_sigs = 0;
            for ant_idx=1:1:num_ants
                if(oversamp_fac>1)
                    curr_samps = resample(rx_interfered_waveforms_babf{ant_idx},1,oversamp_fac);
                else
                    curr_samps = rx_interfered_waveforms_babf{ant_idx};
                end
                summed_sigs = summed_sigs + curr_samps*switch_mat(usr_idx,rfc_idx,ant_idx);
                
            end
            
            summed_sigs = exp(1j*phase_vec(usr_idx,rfc_idx))*summed_sigs;
            [num_samps,~] = size(summed_sigs);
            code_vec = zeros([num_samps,1]);
            code_vec(num_rfc*(usr_idx-1)+rfc_idx:num_rfc*num_users:end) = 1;
            coded_sig = summed_sigs.*code_vec;
            undersamp_samples = summed_sigs(num_rfc*(usr_idx-1)+rfc_idx:num_rfc*num_users:end);
            delay_compensated_samps = delayseq2(undersamp_samples,(num_rfc*(usr_idx-1)+rfc_idx-1)/(num_rfc*num_users));
            
            
            figure(88)
            freq_vec = linspace(-num_rfc*num_users*0.5,num_rfc*num_users*0.5,num_samps);
            subplot(2,1,1)
            abs_vec = abs(fftshift(fft(summed_sigs)))/sqrt(num_samps);
            plot(freq_vec,abs_vec)
            subplot(2,1,2)
            mask_vec = abs_vec>0.1*max(abs_vec);
            plot(freq_vec,angle(fftshift(fft(summed_sigs))).*mask_vec)

            figure(89)
            freq_vec = linspace(-num_rfc*num_users*0.5,num_rfc*num_users*0.5,num_samps);
            subplot(2,1,1)
            abs_vec = abs(fftshift(fft(coded_sig)))/sqrt(num_samps);
            plot(freq_vec,abs_vec)
            subplot(2,1,2)
            mask_vec = abs_vec>0.1*max(abs_vec);
            plot(freq_vec,angle(fftshift(fft(coded_sig))).*mask_vec)
            
            figure(90)
            [num_samps,~] = size(undersamp_samples);
            freq_vec = linspace(-0.5,0.5,num_samps);
            subplot(2,1,1)
            abs_vec = abs(fftshift(fft(undersamp_samples)))/sqrt(num_samps);
            plot(freq_vec,abs_vec)
            subplot(2,1,2)
            mask_vec = abs_vec>0.1*max(abs_vec);
            plot(freq_vec,angle(fftshift(fft(undersamp_samples))).*mask_vec)
            
            figure(91)
            [num_samps,~] = size(undersamp_samples);
            freq_vec = linspace(-0.5,0.5,num_samps);
            subplot(2,1,1)
            abs_vec = abs(fftshift(fft(delay_compensated_samps)))/sqrt(num_samps);
            plot(freq_vec,abs_vec)
            subplot(2,1,2)
            mask_vec = abs_vec>0.1*max(abs_vec);
            plot(freq_vec,angle(fftshift(fft(delay_compensated_samps))).*mask_vec)
            title("Debug 1")
            
            phase_combined_sigs = phase_combined_sigs + delay_compensated_samps;
        end
        kmeans_rate_samples{usr_idx} = phase_combined_sigs;
    end
    
    %  
    
%     for usr_idx=1:1:num_users
%         phase_combined_sigs = 0;
%         for rfc_idx=1:1:num_rfc
%             summed_sigs = 0;
%             for ant_idx=1:1:num_ants
%                 curr_samps = resample(rx_interfered_waveforms_babf{ant_idx},1,meta_params.oversamp_fac);
%                 summed_sigs = summed_sigs + curr_samps*switch_mat(usr_idx,rfc_idx,ant_idx);
%             end
% 
% %             undersamp_samples = exp(-1j*phase_vec(usr_idx,rfc_idx))*summed_sigs(num_rfc*(usr_idx-1)+rfc_idx:num_rfc*num_users:end);
% %             delay_compensated_samps = delayseq2(undersamp_samples,(num_rfc*(usr_idx-1)+rfc_idx-1)/(num_rfc*num_users));
%             phase_combined_sigs = phase_combined_sigs + summed_sigs*exp(1j*phase_vec(usr_idx,rfc_idx));
%         end
%         kmeans_rate_samples{usr_idx} = phase_combined_sigs;
%     end

end


function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins.'));
end
