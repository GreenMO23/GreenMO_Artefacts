function [eval_params,ret_val] = cancel_interf_equalize(meta_params,summed_sigs, ofdm_tx_structs,channel_est_struct)
    
    plot_debug = meta_params.plot_setting;
%     plot_debug_lts_verif = meta_params.plot_setting;
    plot_debug_lts_verif = false;
    
    ofdm_params = ofdm_tx_structs{1};
    ofdm_params.DO_DECODE = 1;
    ofdm_params.DO_APPLY_CFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_SFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_PHASE_ERR_CORRECTION = 1;
    num_users = meta_params.num_users;
    num_lts_to_use = 10;
    [~,num_rfc] = size(summed_sigs);
    lts_thresh = 0.6;
    lts_predetec = true;
    clear op_structs_cell ofdm_params_out
    %% Synch and find lts
    for rfc_idx=1:1:num_rfc
        if(~lts_predetec)
            if(rfc_idx==1)
                [op_struct_cell{rfc_idx}, ofdm_params_out{rfc_idx}] = ofdm_rx(summed_sigs{rfc_idx}, ofdm_params,[], lts_thresh);
                if(isempty(op_struct_cell{rfc_idx}))
                    disp("LTS Corr not found: cancel and equalize")
                    eval_params = -1;
                    ret_val=-1;
                    return
                else
                    ret_val=0;
                    lts_inds = op_struct_cell{1}.lts_ind;
    %                 lts_inds = lts_inds +1; 
                end
            else
                [op_struct_cell{rfc_idx}, ofdm_params_out{rfc_idx}] = ofdm_rx(summed_sigs{rfc_idx}, ofdm_params,lts_inds);
            end
        else
            ret_val =0;
            [op_struct_cell{rfc_idx}, ofdm_params_out{rfc_idx}] = ofdm_rx(summed_sigs{rfc_idx}, ofdm_params,channel_est_struct.lts_inds);
        end
    end
    
    clear op_struct
    
    for rfc_idx=1:1:num_rfc
        op_struct(rfc_idx) =op_struct_cell{rfc_idx};
    end
    %% Extract channel
    conj_idx=1;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct(conj_idx).rx_H_est_vec);
    if(num_packets==1)
        disp("Only 1 packet detected; returning")
        eval_params = -1;
        ret_val = -1;
        return;
    else
        ret_val=0;
    end
    
    ref_chan = op_struct(conj_idx).rx_H_est_vec(:,:,:,start_offset:end);
    ref_mult = exp(-1j*angle(ref_chan));
    
    conj_avg_chans = zeros(num_rfc,num_users,num_subc,num_packets);
    for rfc_idx=1:1:num_rfc
        curr_chan = op_struct(rfc_idx).rx_H_est_vec(:,:,:,start_offset:end);
        conj_chan = squeeze(curr_chan.*ref_mult);
        conj_chan_lts_avg = squeeze(mean(conj_chan(:,ofdm_params.NUM_LTS-num_lts_to_use+1:ofdm_params.NUM_LTS,:,:),2));
        conj_avg_chans(rfc_idx,:,:,:) = conj_chan_lts_avg;
    end
    
    
    %% Segregate LTS, payload for easy time afterwards
    
    sing_lts_len = num_subc; 
    
    lts_len_singuser = (((ofdm_params.NUM_LTS+0.5)*sing_lts_len)); % 10 LTS + 0.5 CP
    lts_len = (((ofdm_params.NUM_LTS+0.5)*sing_lts_len)+ofdm_params.inter_user_zeros)*num_users; % total preamb len
    payload_len = (ofdm_params.N_SC+ofdm_params.CP_LEN)*ofdm_params.N_OFDM_SYMS+0.5*ofdm_params.inter_user_zeros;
    payload_len_ceiled = ceil(payload_len/(ofdm_params.N_SC+ofdm_params.CP_LEN))*(ofdm_params.N_SC+ofdm_params.CP_LEN);
    pkt_and_payload_plus_zeros = payload_len_ceiled+lts_len; % +0.5*ofdm_params.inter_user_zeros for multipath
    num_subc_groups = pkt_and_payload_plus_zeros/ofdm_params.N_SC;
    
    % RX_LTS stores all the recovered LTS
    rx_lts = zeros(num_packets,ofdm_params.NUM_LTS,num_users,num_rfc,num_subc);
    % Payload and LTSes for SISO processing
    rx_payload_only = zeros(num_packets, num_rfc,payload_len_ceiled);
    rx_all_ltses = zeros(num_packets, num_rfc,lts_len);
    rx_payload_and_lts = zeros(num_packets, num_rfc,pkt_and_payload_plus_zeros);
    
    for pkt_idx=1:1:num_packets
       for rfc_idx=1:1:num_rfc
           lts_offsets = op_struct(rfc_idx).user_lts_offsets;
           pkt_lts_ind = op_struct(rfc_idx).packet_lts_list(pkt_idx);
           pkt_payload_ind = op_struct(rfc_idx).packet_payload_list(pkt_idx);
           raw_samps_rfi = summed_sigs{rfc_idx};
           rx_payload_and_lts(pkt_idx,rfc_idx,:) = raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
           rx_payload_only(pkt_idx,rfc_idx,:) = raw_samps_rfi(pkt_lts_ind+lts_len:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
           rx_all_ltses(pkt_idx,rfc_idx,:) = raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+lts_len-1);
%            if(plot_debug)
%                subplot(floor(num_users/2),2,usr_idx)
%                xline(pkt_lts_ind,'r--');
%                xline(pkt_payload_ind,'g-.');
%                xline(pkt_lts_ind+pkt_and_payload_plus_zeros,'m:');
%            end
           all_lts = raw_samps_rfi(pkt_lts_ind:pkt_payload_ind-1);
           rel_lts_ind = op_struct(rfc_idx).packet_lts_list(pkt_idx,:)-pkt_lts_ind+1;
           lts_offset=0;
           half_lts = num_subc/2;
           if(plot_debug_lts_verif)
               figure(10)
               plot(abs(all_lts))
               hold on
           end
           lts_start_ind = 1;
           for user_idx=1:1:num_users
               lts_start_ind_offseted = lts_start_ind+lts_offsets(user_idx);
               for lts_idx=1:1:ofdm_params.NUM_LTS
                   if(plot_debug_lts_verif)
                      xline(lts_start_ind_offseted+half_lts+(lts_idx-1)*num_subc,'r:');
                   end
    %                 (num_packets,ofdm_params.NUM_LTS,num_users,num_rfc,num_subc)
                   rx_lts(pkt_idx,lts_idx,user_idx,rfc_idx,:) = all_lts(lts_start_ind_offseted+half_lts+(lts_idx-1)*num_subc:lts_start_ind_offseted+half_lts+(lts_idx)*num_subc-1);
               end
               
%                lts_start_ind = lts_start_ind+(ofdm_params.NUM_LTS+0.5)*num_subc+ofdm_params.inter_user_zeros;
    %            xline(lts_start_ind,'m--')
           end
       end
    end
    
    %% Intermediate data variables
    pkt_to_use = 1;
    chan_to_use = conj_avg_chans(:,:,:,pkt_to_use);
%     chan_to_use = digital_channel_mat(:,:,:,pkt_to_use);
%         chan_to_use = squeeze(mean(conj_avg_chans,[1]));

    nonnull_subc = [2:27 39:64];
    num_nonnull_subc = numel(nonnull_subc);
    % rx_lts = zeros(num_packets,ofdm_params.NUM_LTS,num_users,num_rfc,num_subc);
    combined_lts_f = zeros(num_packets,num_users,num_users,ofdm_params.NUM_LTS,num_subc);
    user_comb_vecs = zeros(num_users,num_rfc,num_subc);
    tdma_snrs = zeros(num_packets,num_users,num_rfc,num_nonnull_subc);
    % Take FFTCombine and take IFFT to get back to td
%     rx_payload_and_lts = zeros(num_packets, num_rfc,pkt_and_payload_plus_zeros);
%     rx_payload_and_lts_groups_64 = reshape(rx_payload_and_lts,num_packets,num_rfc,num_subc,num_subc_groups);
%     if(plot_debug)
%         figure(81)
%         plot(abs(squeeze(rx_payload_and_lts(1,1,:))));
%         hold on
%         all_rx_lts=[];
%         for i=1:1:num_users
%             all_rx_lts = [all_rx_lts  zeros(1,0.5*sing_lts_len)];
%             for j=1:1:ofdm_params.NUM_LTS
%                 % rx_lts = zeros(num_packets,ofdm_params.NUM_LTS,num_users,num_rfc,num_subc);
%                 all_rx_lts = [all_rx_lts squeeze(rx_lts(1,j,i,1,:)).'];
%             end
%             all_rx_lts = [all_rx_lts  zeros(1,ofdm_params.inter_user_zeros)];
%         end
%         all_rx_lts = all_rx_lts .'; 
%         plot([abs(all_rx_lts).' abs(squeeze(rx_payload_only(1,1,:))).']);
%         plot([abs(squeeze(rx_all_ltses(1,1,:))).' abs(squeeze(rx_payload_only(1,1,:))).']);
% %         num_subc_groups
% %         for i=1:1:num_subc_groups
% %             xline(64*i);
% %         end
%     end
    
    num_syms_group = ceil(payload_len/(ofdm_params.N_SC+ofdm_params.CP_LEN));
    rx_payload_only_cp_reshaped = reshape(rx_payload_only,num_packets,num_rfc,num_subc+ofdm_params.CP_LEN,num_syms_group);
    rx_payload_only_cp_rem = rx_payload_only_cp_reshaped(:,:,ofdm_params.CP_LEN+1:num_subc+ofdm_params.CP_LEN,:);
    rx_payload_f = fft(rx_payload_only_cp_rem,64,3); 
    combined_payload_and_lts_f = zeros(num_packets,num_users,num_subc,num_syms_group);
    
    %% Combining codes now finally
    for pkt_idx=1:1:num_packets
        % Get freq domain LTS's for each packet
        if(num_users>1)
            rx_lts_f = fft(squeeze(rx_lts(pkt_idx,:,:,:,:)),num_subc,4);% 4th axis is subcarrier
        else
            rx_lts_f = fft(squeeze(rx_lts(pkt_idx,:,:,:,:)),num_subc,3);% 4th axis is subcarrier
        end
        % perform combining per subcarrier
        for subc_idx=1:1:num_nonnull_subc 
            % Finally obtain narrowband channel matrix for each subc
            if(num_users==1)
                chan_mat = chan_to_use(:,nonnull_subc(subc_idx));
            else
                chan_mat = chan_to_use(:,:,nonnull_subc(subc_idx));
            end
            % Now go on and implement digital combining
            for user_idx=1:1:num_users
                if(num_users>1)
                    %% MMSE postcoding
                    noise_rat = 10^(30/10);
                    H_hat = chan_mat;
%                     norm_H_hat=H_hat/norm(H_hat,'fro');
%                     inv_mat=inv(norm_H_hat*(norm_H_hat')+eye(num_users)/noise_rat);
%                     W=(norm_H_hat.')*inv_mat;
                    if(meta_params.skip_combining)
                        W=eye(num_users);
                    else
                        W = pinv(H_hat);
                    end
%                     W=W/norm(W,'fro');
                    fin_comb_vec = W(user_idx,:).';
                    user_comb_vecs(user_idx,:,nonnull_subc(subc_idx)) = fin_comb_vec;
                    %% Combine with W
                    
                    curr_payload_to_be_combined = squeeze(rx_payload_f(pkt_idx,:,nonnull_subc(subc_idx),:)).';
                    repmat_combvec = repmat(fin_comb_vec.',[num_syms_group,1]);
                    combined_samp = sum(curr_payload_to_be_combined.*repmat_combvec,2); % RFCs are getting summed over here
                    combined_payload_and_lts_f(pkt_idx,user_idx,nonnull_subc(subc_idx),:)=combined_samp;
                    
                    %% Compute TDMA SNRs from uncombined LTS's in freq domain (For each packet and each subcarrier)
                    curr_lts_tdmasnr = squeeze(rx_lts_f(:,user_idx,:,nonnull_subc(subc_idx)));
                    sig_pow = squeeze(mean(abs(curr_lts_tdmasnr).^2,1));
                    noise_pow = squeeze(var(abs(curr_lts_tdmasnr),0,1));
                    tdma_snrs(pkt_idx,user_idx,:,subc_idx) = sig_pow./noise_pow;

                    for comb_user_idx=1:1:num_users
                        curr_lts_to_be_combined = squeeze(rx_lts_f(:,comb_user_idx,:,nonnull_subc(subc_idx)));
                        repmat_combvec = repmat(fin_comb_vec.',[ofdm_params.NUM_LTS,1]);
                        combined_samp = sum(curr_lts_to_be_combined.*repmat_combvec,2);
                        combined_lts_f(pkt_idx,user_idx,comb_user_idx,:,nonnull_subc(subc_idx))=combined_samp;     
                    end
                else
                    %% Just populat stuff for single user case, no need for any coding
                    ref_chain_idx =1;
                    curr_payload_to_be_combined = squeeze(rx_payload_and_lts_groups_64_f(pkt_idx,ref_chain_idx,nonnull_subc(subc_idx),:)).';
                    curr_lts_to_be_combined = squeeze(rx_lts_f(:,ref_chain_idx,nonnull_subc(subc_idx)));

                    sig_pow = squeeze(mean(abs(curr_lts_to_be_combined).^2,1));
                    noise_pow = squeeze(var(abs(curr_lts_to_be_combined),0,1));
                    tdma_snrs(pkt_idx,user_idx,:,subc_idx) = sig_pow./noise_pow;

                    combined_lts_f(pkt_idx,user_idx,user_idx,:,nonnull_subc(subc_idx))=curr_lts_to_be_combined;
                    combined_payload_and_lts_f(pkt_idx,user_idx,nonnull_subc(subc_idx),:)=curr_payload_to_be_combined;
                end
            end
        end
    end
    
    %% Analysis, EVM BER and all those shizz
    % get back to time domain and flatten the array
    combined_lts_t = ifft(combined_lts_f,num_subc,5);
    combined_payload_t = ifft(combined_payload_and_lts_f,num_subc,3);
    combined_payload_t_plus_cp = zeros(num_packets,num_users,num_subc+ofdm_params.CP_LEN,num_syms_group);
    combined_payload_t_plus_cp(:,:,1:ofdm_params.CP_LEN,:) = combined_payload_t(:,:,num_subc-ofdm_params.CP_LEN+1:num_subc,:);
    combined_payload_t_plus_cp(:,:,ofdm_params.CP_LEN+1:ofdm_params.CP_LEN+num_subc,:) = combined_payload_t;
    combined_payload_fin = reshape(combined_payload_t_plus_cp,num_packets,num_users,payload_len_ceiled,[]);
%     plot_debug=true;
    if(plot_debug)
        pkt_idx = 1;
        figure(83)
        num_rows=2;
        num_cols=ceil(num_users/num_rows); 
        for i=1:1:num_users
            subplot(num_cols,num_rows,i)
            all_rx_lts=[];
            for ii=1:1:num_users
                all_rx_lts = [all_rx_lts  zeros(1,0.5*sing_lts_len)];
                for j=1:1:ofdm_params.NUM_LTS
                    % combined_lts_f = zeros(num_packets,num_users,num_users,ofdm_params.NUM_LTS,num_subc);
                    % combined_lts_f(pkt_idx,user_idx,comb_user_idx,:,nonnull_subc(subc_idx))=combined_samp;     
                    all_rx_lts = [all_rx_lts squeeze(combined_lts_t(pkt_idx,i,ii,j,:)).'];
                end
                all_rx_lts = [all_rx_lts  zeros(1,ofdm_params.inter_user_zeros)];
            end
            all_rx_lts = all_rx_lts .'; 
            
            
%             plot([abs(squeeze(rx_all_ltses(pkt_idx,i,:))).' abs(squeeze(rx_payload_only(pkt_idx,i,:))).']);
%             hold on
            plot([abs(all_rx_lts).' squeeze(abs(combined_payload_fin(pkt_idx,i,:))).']);
        end
    end
                                             % pkt,user,user,10lts,nonnulls
    combined_lts_f_non_null = combined_lts_f(:,:,:,:,nonnull_subc);
    
    % these are the things you need
    sinr_calc = zeros(num_packets,num_users,num_nonnull_subc );
    ber_calc = zeros(num_packets,num_users);
    evm_snr_calc = zeros(num_packets,num_users);
    
    clear all_users_decode_struct
    for pkt_idx=1:1:num_packets
        clear decode_params decode_op_structs
        for user_idx=1:1:num_users
            subc_pow_profile = (abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:))).^2);
            if(num_users>1)
                mean_pow_per_comb_user = squeeze(mean(subc_pow_profile,2)); %mean across 10lts
                subc_mag_profile = abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:)));
                noise_pow = squeeze(var(subc_mag_profile(user_idx,:,:),0,2));
                sig_pow = mean_pow_per_comb_user(user_idx,:).';
                curr_user = user_idx;
                user_vec = 1:1:num_users;
                interf_users = user_vec;
                interf_users(curr_user) =[];
                interf_pow = sum(mean_pow_per_comb_user(interf_users,:),1).';
                if(meta_params.skip_combining)
                	sinr_vec = sig_pow./(noise_pow);
                else
                    sinr_vec = sig_pow./(interf_pow+noise_pow);
                end
                sinr_calc(pkt_idx,user_idx,:) = sinr_vec;
            else
                mean_pow_per_comb_user = squeeze(mean(subc_pow_profile,1)); %mean across 10lts
                subc_mag_profile = abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:)));
                noise_pow = squeeze(var(subc_mag_profile,0,1));rx_waveforms_per_user, 
                sig_pow = mean_pow_per_comb_user(user_idx,:);
                sinr_vec = sig_pow./(noise_pow);
                sinr_calc(pkt_idx,user_idx,:) = sinr_vec;
            end    
            
            
            all_rx_lts=[];
            for ii=1:1:num_users
                all_rx_lts = [all_rx_lts  zeros(1,0.5*sing_lts_len)];
                for j=1:1:ofdm_params.NUM_LTS
                    % combined_lts_f = zeros(num_packets,num_users,num_users,ofdm_params.NUM_LTS,num_subc);
                    % combined_lts_f(pkt_idx,user_idx,comb_user_idx,:,nonnull_subc(subc_idx))=combined_samp;     
                    all_rx_lts = [all_rx_lts squeeze(combined_lts_t(pkt_idx,user_idx,ii,j,:)).'];
                end
                all_rx_lts = [all_rx_lts  zeros(1,ofdm_params.inter_user_zeros)];
            end
            all_rx_lts = all_rx_lts .'; 
            curr_payload = [all_rx_lts.' squeeze(combined_payload_fin(pkt_idx,user_idx,:)).'];
            ofdm_params_tx = ofdm_tx_structs{user_idx};
            [decode_op_structs(user_idx),decode_params(user_idx)]= ofdm_decode(squeeze(curr_payload).', ofdm_params_tx, ofdm_params,user_idx-1, lts_offsets);
%             [decode_op_structs(user_idx),decode_params(user_idx)]= ofdm_decode(summed_sigs{user_idx}, ofdm_params_tx, ofdm_params,user_idx-1, lts_offsets);
            ber_calc(pkt_idx, user_idx) = decode_op_structs(user_idx).berratio;
            evm_snr_calc(pkt_idx, user_idx) = decode_op_structs(user_idx).snr;
        end
        all_users_decode_struct{pkt_idx} = decode_op_structs;
    end
    
    
    eval_params.ber_calc = ber_calc;
    eval_params.evm_snr_calc = 10*log10(evm_snr_calc);
    eval_params.sinr_calc = 10*log10(squeeze(mean(sinr_calc,[3])));
    eval_params.tdma_snrs = 10*log10(squeeze(mean(tdma_snrs,[1,4])));
    
    
%     evm_snr_packets = mean(evm_snr_calc,2).';
%     mean_snr = mean(evm_snr_packets);
%     std_snr =  std(evm_snr_packets);
%     tol = 2;
%     good_snrs = (abs(evm_snr_packets-mean_snr)<tol*std_snr);
%     cleaned_snrs = evm_snr_packets(find(good_snrs));
%     avg_evm_snr_across_users = [avg_evm_snr_across_users cleaned_snrs];
% 
%     avg_evm_snr_users_across_packets(test_index,:) = 10*log10(mean(evm_snr_calc,1));
%     min_evm_snr_users_across_packets(test_index,:) = 10*log10(min(evm_snr_calc,[],1));
%     
%     sinr_subc_averaged = squeeze(mean(sinr_calc,[3])); 
%     avg_sinr_users_across_packets(test_index,:) = 10*log10(mean(sinr_subc_averaged,1));
%     min_sinr_users_across_packets(test_index,:) = 10*log10(min(sinr_subc_averaged,[],1));
%     
%     avg_ber_users_across_packets(test_index,:) = mean(ber_calc,1);
%     max_ber_users_across_packets(test_index,:) = min(ber_calc,[],1);
    
end

