function [chan_est_struct,snr_db,ret_val] = synch_and_get_channels_twouser(meta_params,ofdm_tx_structs,channel_struct,rx_interfered_waveforms)
    
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
    if(downsample_flag==0)
        downsamp_conj = resample(rx_interfered_waveforms{conj_ant},1,oversamp_fac);
    else
        downsamp_conj = downsamp_filt(rx_interfered_waveforms{conj_ant},oversamp_fac);
    end
    [op_struct_conj, ~] = ofdm_rx(downsamp_conj, ofdm_params_to_use);

    if(isempty(op_struct_conj))
        disp("LTS not detected, exiting")
        chan_est = -1;
        ret_val =-1;
        chan_est_struct=-1; 
        snr_db=-1;
        return;
    else
        ret_val =0;
    end
    
    lts_inds = op_struct_conj.lts_ind;
    chan_est_struct.lts_inds = lts_inds;

    clear downsampled_sigs 
    for ant_ind=1:1:num_ants
        if(downsample_flag==0)
            downsamp_samp = resample(rx_interfered_waveforms{ant_ind},1,oversamp_fac);
        else
            downsamp_samp = downsamp_filt(rx_interfered_waveforms{ant_ind},oversamp_fac);
        end
        downsampled_sigs{ant_ind} = downsamp_samp;
        % call ofdm rx with lts inds
        [op_struct_cell{ant_ind}, ~] = ofdm_rx(downsamp_samp, ofdm_params_to_use, lts_inds);
    end

    op_struct = op_struct_cell;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct{1}.rx_H_est_vec);
    if(num_packets==1)
        disp("Only 1 packet detected; returning")
        chan_est = -1;
        ret_val =-1;
        chan_est_struct=-1; 
        snr_db=-1;
        return;
    else
        ret_val=0;
    end
    conj_avg_chans = zeros(num_packets,num_users,num_ants,num_subc);
%     conj_idx =3;
    ref_chan = op_struct{conj_ant}.rx_H_est_vec(1,:,:,start_offset:end);
    ref_mult = exp(-1j*angle(ref_chan));

    partial_ant_vec=1:1:64;
    partial_ant_vec = reshape(partial_ant_vec,[8,8]).';
%     chan_est_struct.hyb_phase_angles = zeros([1,num_ants]);
    
    for user_idx=1:1:num_users
        if(disp_plots)
            figure(user_idx+fig_idx)
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
%                 
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
    
    
    chan_est_struct.num_lts_to_use = num_lts_to_use;
    chan_est_struct.num_ants = num_ants;
    

    babf_chan_mat = squeeze(conj_avg_chans(1,:,:,4));
    if(num_ants == oversamp_fac)
        chan_est_struct.babf_switched_mat = eye(num_ants);
    else
        babf_switched_mat = analog_beamf_config_choser(babf_chan_mat',num_ants,1:1:num_ants,pi/2);
        chan_est_struct.babf_switched_mat = babf_switched_mat;
    end

    
%     chan_est_struct.babf_switched_mat = babf_switched_mat;
%     
%     filledln_switched_mat = filledln_config_choser(babf_chan_mat');
%     chan_est_struct.filledln_switched_mat = filledln_switched_mat;
    
%     partial_ant_vec=1:1:64;
%     partial_ant_vec = reshape(partial_ant_vec,[8,8]).';
    
%     babf_chan_mat_norm = babf_chan_mat./abs(babf_chan_mat);
    chan_est_struct.fully_connected_mat = babf_chan_mat';
%     chan_est_struct.partially_connected_mat = babf_chan_mat';
    
    
%     nonnull_subc = [2:27 39:64];
%     kmeans_chan_mat = squeeze(conj_avg_chans(1,:,:,nonnull_subc));
%     chan_est_struct.kmeans_chan_mat = squeeze(conj_avg_chans(1,:,:,nonnull_subc));
%     chan_est_struct.kmeans_chan_mat = squeeze(conj_avg_chans(1,:,:,[2,3]));
%      kmeans_chan_mat = squeeze(conj_avg_chans(1,:,:,[2,3]));
%     chan_est_struct.kmeans_switched_mat = kmeans_config_choser(kmeans_chan_mat,1);
end

function kmeans_switched_mat = kmeans_config_choser(kmeans_chan_mat,num_slots)
    [num_users,num_ants,num_subc]=size(kmeans_chan_mat);
    num_classes = 3;
    if(num_slots==1)
        switched_mat = zeros(num_ants,num_users);
        kmeans_class_vec = zeros(num_ants,num_users);
        dominant_class = zeros(1,num_users);
        class_counts = zeros(num_users,num_classes);
        for user_idx=1:1:num_users
            curr_chan_mat = squeeze(kmeans_chan_mat(user_idx,:,:));
            kmeans_class = kmeans(angle(curr_chan_mat),num_classes).';
            kmeans_class_vec(:,user_idx) = kmeans_class;
            dominant_class(user_idx) = mode(kmeans_class);
            switched_mat(:,user_idx) = kmeans_class==mode(kmeans_class);
        end
        rank_val = rank(switched_mat);
        if(rank_val==num_users) % exit out of the loop if rank is full
            kmeans_switched_mat = switched_mat;
            return
        else
            disp("Finding perturbed k means solns")
            non_invert_flag = true;
            curr_perturb_user_idx = 2;
            global_perturb_idx=1;
            class_vec = dominant_class;
            while(non_invert_flag)
                class_to_select = mod(dominant_class(curr_perturb_user_idx)+1,num_classes)+1;
                class_vec(curr_perturb_user_idx) = class_to_select;
                switched_mat(:,curr_perturb_user_idx) = kmeans_class==class_to_select;
                rank_val = rank(switched_mat);
                if(rank_val==num_users)
                    kmeans_switched_mat = switched_mat;
                    return
                else
                    disp("Finding perturbed k means solns: lvl 2")
                    class_vec(curr_perturb_user_idx) = dominant_class(curr_perturb_user_idx);
                    curr_perturb_user_idx = mod(curr_perturb_user_idx +1,num_users)+1;
                    if(curr_perturb_user_idx==2)
                        dominant_class(global_perturb_idx) = mod(dominant_class(global_perturb_idx)+1,num_classes);
                        global_perturb_idx= mod(curr_perturb_user_idx +1,num_users)+1;
                    end
                    
                    if(global_perturb_idx==num_users)
                        % choose random matrix
                        rand_flag = true;
                        while(rand_flag)
                            switched_mat= randi(2,num_ants,num_users)-1; 
                            rank_val = rank(switched_mat);
                            if(rank_val==num_users)
                                kmeans_switched_mat = switched_mat;
                                return
                            end
                        end
                    end
                end    
            end
        end
    else
        kmeans_switched_mat = -1;
        return
    end
end

function switched_beamf_mat = analog_beamf_config_choser(beamf_mat_all_users,num_ants_to_use,ants_to_use,inphase_val)
    [max_ants,num_users] = size(beamf_mat_all_users);
    mat_to_use = beamf_mat_all_users(ants_to_use,:);
    non_invert_chan_mat = true;
    ant_idx_perf=1;
    perf_choice = ones(1,num_users); % choose the best ones if possible
    while(non_invert_chan_mat)
        switched_beamf_mat= zeros(max_ants,num_users);
        if(sum(perf_choice>ceil(max_ants/2)+1)>0)
            break
        end
        for user_idx=1:1:num_users
            antenna_score_func = zeros(1,num_ants_to_use);
            antenna_configs = zeros(num_ants_to_use,num_ants_to_use);
            curr_vec = mat_to_use(:,user_idx); % user's channel 
            for ant_idx=1:1:num_ants_to_use % iterate over all antennas
                curr_phase = mat_to_use(ant_idx,user_idx);
                conj_vec = curr_vec*conj(curr_phase); % conjugate multiplication allows to compute angle
                antenna_configs(ant_idx,:) = abs(angle(conj_vec))<inphase_val; % see if angle is less than the desired inphase angle
                antenna_score_func(ant_idx) = -sum(antenna_configs(ant_idx,:)); % sort does asc, want to choose most number of antennas
            end
            [~,antenna_perf_order] = sort(antenna_score_func); % sort such that my configs which show highest numnber of antennas are on top
            switched_beamf_mat(ants_to_use,user_idx) = antenna_configs(antenna_perf_order(perf_choice(user_idx)),:); % choose from the perf order
        end
        
        rank_val = rank(switched_beamf_mat);
        if(rank_val==num_users) % exit out of the loop if rank is full
            non_invert_chan_mat=false;
            break
        else
            perf_choice(ant_idx_perf)=perf_choice(ant_idx_perf)+1;
            ant_idx_perf=mod((ant_idx_perf+1)-1,num_users)+1;
        end
            
    end
    
    
    while(non_invert_chan_mat)
%         if(num_ants_to_use==4)
%             switched_beamf_mat= zeros(max_ants,num_users);
%             switched_beamf_mat(ants_to_use,:) = fliplr(eye(4));
%             disp("Unable to find invertible analog soln, supplying identity mat")
%             break
%         end
        switched_beamf_mat= zeros(max_ants,num_users);
        switched_beamf_mat(ants_to_use,:) = randi(2,num_ants_to_use,num_users)-1;
        if(rank(switched_beamf_mat)==num_users)
            non_invert_chan_mat=false;
            disp("Unable to find invertible analog soln, supplying random mat")
            break
        end
    end

    
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

