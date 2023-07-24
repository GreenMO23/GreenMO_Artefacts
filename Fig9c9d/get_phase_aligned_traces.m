function [rx_lts,rx_payload_only,rx_all_ltses,rx_payload_and_lts, chan_to_use,lts_offsets,ret_val] = get_phase_aligned_traces(rx_samps_ch_est,ofdm_params)

%     traces{1} = rx_samps_ch_est{2};
%     ant1345_traces = rx_samps_ch_est{3};
%     ant1678_traces = rx_samps_ch_est{4};
    
    rfc_index = 2; % corresponding to warp
%     conj_index = 4;
    
    num_users = 4;
    samp_fac = 4;
    ofdm_params.DO_DECODE = 1;
    ofdm_params.DO_APPLY_CFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_SFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_PHASE_ERR_CORRECTION = 1;
    num_packets_arr = zeros([1,3]);
    total_rfc = 8;
%     chan_to_use = zeros(num_users,total_rfc,ofdm_params.N_SC);
    chan_to_use = zeros(total_rfc,num_users,ofdm_params.N_SC);
    
    % 3 times data is collected overpkt_to_use = 1;
    clear rx_samples_resamped op_struct_cell    
    for i=1:1:3
        rx_samps = rx_samps_ch_est{1+i};
        rx_samps_rfc = rx_samps(rfc_index,:).';
        
        for oversamp_ind=1:1:samp_fac
            rx_samples_resamped{i,oversamp_ind} = delayseq2(rx_samps_rfc(oversamp_ind:samp_fac:end),(oversamp_ind-1)/samp_fac);
        end
    
        % Obtain LTS indexes and payloads and segregate
        num_rfc = 4;
        rfc_map = [2,3,4,1];
        for rfc_idx=rfc_map
            if(rfc_idx==2)
                [op_struct_cell{i,rfc_idx}, ofdm_params_out] = ofdm_rx(rx_samples_resamped{i,rfc_idx}, ofdm_params);
                if(isempty(op_struct_cell{i,rfc_idx}))
%                     disp("LTS Corr not found, exiting with zeros")
                    ret_val=-1;
                    rx_lts = -1;
                    rx_payload_only=-1;
                    rx_all_ltses=-1;
                    rx_payload_and_lts=-1;
                    chan_to_use=-1;
                    lts_offsets=-1;
                    return
                else
                    ret_val=0;
                    lts_inds = op_struct_cell{i,2}.lts_ind;
    %                 lts_inds = lts_inds +1; 
                end
            else
                [op_struct_cell{i,rfc_idx}, ofdm_params_out] = ofdm_rx(rx_samples_resamped{i,rfc_idx}, ofdm_params,lts_inds);
            end
        end
    
        [~,~,num_subc,num_packets] = size(op_struct_cell{i,1}.rx_H_est_vec);
        num_packets_arr(i) = num_packets;
    end
    
    num_min_packets = min(num_packets_arr);
    
    sing_lts_len = ofdm_params.N_SC; 
    lts_len_singuser = (((ofdm_params.NUM_LTS+0.5)*sing_lts_len)); % 10 LTS + 0.5 CP
    lts_len = (((ofdm_params.NUM_LTS+0.5)*sing_lts_len)+ofdm_params.inter_user_zeros)*num_users; % total preamb len
    payload_len = (ofdm_params.N_SC+ofdm_params.CP_LEN)*ofdm_params.N_OFDM_SYMS+0.5*ofdm_params.inter_user_zeros;
    payload_len_ceiled = ceil(payload_len/(ofdm_params.N_SC+ofdm_params.CP_LEN))*(ofdm_params.N_SC+ofdm_params.CP_LEN);
    pkt_and_payload_plus_zeros = payload_len_ceiled+lts_len; % +0.5*ofdm_params.inter_user_zeros for multipath
    num_subc_groups = pkt_and_payload_plus_zeros/ofdm_params.N_SC;
    
    
    
    %% important stuff
    % RX_LTS stores all the recovered LTS, per user
    rx_lts = zeros(num_min_packets,ofdm_params.NUM_LTS,num_users,total_rfc,num_subc);
    debug_lts = zeros(num_min_packets,ofdm_params.NUM_LTS,num_users,total_rfc,num_subc);
    % Only the payloads
    rx_payload_only = zeros(num_min_packets, total_rfc,payload_len_ceiled);
    % RX_LTS stores all the LTSes
    rx_all_ltses = zeros(num_min_packets, total_rfc,lts_len);
    rx_all_ltses_comp = zeros(num_min_packets, total_rfc,lts_len);
    % Payload and LTSes for SISO processing
    rx_payload_and_lts = zeros(num_min_packets, total_rfc,pkt_and_payload_plus_zeros);
    
%     rx_lts = zeros(num_min_packets,ofdm_params.NUM_LTS,num_users,total_rfc,num_subc);
    % Only the payloads
    ref_rx_payload_only = zeros(num_min_packets,3, payload_len_ceiled);
    % RX_LTS stores all the LTSes
    ref_rx_all_ltses = zeros(num_min_packets,3, lts_len);
    % Payload and LTSes for SISO processing
    ref_rx_payload_and_lts = zeros(num_min_packets,3, pkt_and_payload_plus_zeros);
    
    
    
    %% antenna mapping
    conj_idx = 2;

    ant_chan_idx = zeros(1,total_rfc);
    ant_chan_num_idx = zeros(1,total_rfc);
    ref_ant_idx = conj_idx;
    ant_chan_idx(1) = conj_idx;
    zero_chan_idx_1 = mod((conj_idx+1)-1,4)+1;
    ant_chan_idx(2) = mod((conj_idx+2)-1,4)+1;
    ant_chan_idx(3) = mod((conj_idx+3)-1,4)+1;
    ant_chan_num_idx(1:3) = ones(1,3);
    ant_chan_idx(4) = mod((conj_idx+1)-1,4)+1;
    ant_chan_idx(5) = mod((conj_idx+2)-1,4)+1;
    ant_chan_idx(6) = mod((conj_idx+3)-1,4)+1;
    ant_chan_num_idx(4:6) = 2*ones(1,3);
    ant_chan_idx(7) = mod((conj_idx+1)-1,4)+1;
    ant_chan_idx(8) = mod((conj_idx+2)-1,4)+1;
    zero_chan_idx_2 = mod((conj_idx+3)-1,4)+1;
    ant_chan_num_idx(7:8) = 3*ones(1,2);
    
    ants_to_select = {[1,2,3],[4,5,6],[7,8]};
    num_ants_to_select = [3,3,2];
    num_lts_to_use = 10;
    
    
    for i=1:1:3
        
        clear op_struct

        for rfc_idx=1:1:num_rfc
            op_struct(rfc_idx) =op_struct_cell{i,rfc_idx};
        end
        
        %% channel extraction
        start_offset = 1;
        [~,~,~,num_packets] = size(op_struct(ref_ant_idx).rx_H_est_vec);
        ref_chan = op_struct(ref_ant_idx).rx_H_est_vec(:,:,:,start_offset:end);
        ref_mult = exp(-1j*angle(ref_chan));
        conj_avg_chans = zeros(num_ants_to_select(i),num_users,num_subc,num_packets);
        ant_indexes = ant_chan_idx(ants_to_select{i});
        global_ant_indexes = ants_to_select{i};
        
        for rfc_idx=1:1:num_ants_to_select(i)
            ant_idx_wrapped = ant_indexes(rfc_idx); 
            curr_chan = op_struct(ant_idx_wrapped).rx_H_est_vec(:,:,:,start_offset:end);
            conj_chan = squeeze(curr_chan.*ref_mult);
            conj_chan_lts_avg = squeeze(mean(conj_chan(:,ofdm_params.NUM_LTS-num_lts_to_use+1:ofdm_params.NUM_LTS,:,:),2));
            conj_avg_chans(rfc_idx,:,:,:) = conj_chan_lts_avg;
        end
        
        % take channel for packet 1
        chan_to_use(global_ant_indexes,:,:) = conj_avg_chans(:,:,:,1);
        
        plot_debug_lts_verif = false;
        
        for pkt_idx=1:1:num_min_packets
            ant_indexes = ant_chan_idx(ants_to_select{i});
            global_ant_indexes = ants_to_select{i};
            
            for rfc_idx=1:1:num_ants_to_select(i)
               
                
               ant_idx_wrapped = ant_indexes(rfc_idx); 
               global_ant_idx = global_ant_indexes(rfc_idx);
                
               lts_offsets = op_struct(ant_idx_wrapped).user_lts_offsets;
               pkt_lts_ind = op_struct(ant_idx_wrapped).packet_lts_list(pkt_idx);
               raw_samps_rfi = rx_samples_resamped{i,ant_idx_wrapped};
               
               curr_payload_and_lts = raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
               curr_payload_only = raw_samps_rfi(pkt_lts_ind+lts_len:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
               curr_all_ltses = raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+lts_len-1);
               
%                if(i==1 && rfc_idx==1)
                if(rfc_idx==1)
                    refi_raw_samps_rfi = rx_samples_resamped{i,ref_ant_idx};
                    pkt_lts_ind = op_struct(ref_ant_idx).packet_lts_list(pkt_idx);
                    ref_rx_payload_and_lts(pkt_idx,i,:) = refi_raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
                    ref_rx_payload_only(pkt_idx,i,:) = refi_raw_samps_rfi(pkt_lts_ind+lts_len:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
                    ref_rx_all_ltses(pkt_idx,i,:) = refi_raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+lts_len-1);
                end
%               
%                elseif(rfc_idx==1)
%                    refi_lts_offsets = op_struct(ref_ant_idx).user_lts_offsets;
%                    refi_pkt_lts_ind = op_struct(ref_ant_idx).packet_lts_list(pkt_idx);
%                    refi_raw_samps_rfi = rx_samples_resamped{i,ref_ant_idx};
% 
%                    refi_payload_and_lts = refi_raw_samps_rfi(refi_pkt_lts_ind:refi_pkt_lts_ind+pkt_and_payload_plus_zeros-1);
%                    refi_payload_only = refi_raw_samps_rfi(refi_pkt_lts_ind+lts_len:refi_pkt_lts_ind+pkt_and_payload_plus_zeros-1);
%                    refi_all_ltses = refi_raw_samps_rfi(refi_pkt_lts_ind:refi_pkt_lts_ind+lts_len-1);
%                     
%                    conj_mult_payload_lts = squeeze(conj(ref_rx_payload_and_lts(pkt_idx,1,:))).';
%                    conj_mult_payload_only = squeeze(conj(ref_rx_payload_only(pkt_idx,1,:))).';
%                    conj_mult_all_ltses = squeeze(conj(ref_rx_all_ltses(pkt_idx,1,:))).';
%                    
%                    ref_rx_payload_and_lts(pkt_idx,i,:) = exp(1j*angle(refi_payload_and_lts.*conj_mult_payload_lts));
%                    ref_rx_payload_only(pkt_idx,i,:) = exp(1j*angle(refi_payload_only.*conj_mult_payload_only));
%                    ref_rx_all_ltses(pkt_idx,i,:) = exp(1j*angle(refi_all_ltses.*conj_mult_all_ltses));
%                    
%                end
               
%                if(i~=1)
               conj_mult_payload_lts = squeeze(conj(ref_rx_payload_and_lts(pkt_idx,i,:))).';
               conj_mult_payload_only = squeeze(conj(ref_rx_payload_only(pkt_idx,i,:))).';
               conj_mult_all_ltses = squeeze(conj(ref_rx_all_ltses(pkt_idx,i,:))).';

               phase_compensation_payload_lts = exp(1j*angle(curr_payload_and_lts.*conj_mult_payload_lts ));
               phase_compensation_payload = exp(1j*angle(curr_payload_only.*conj_mult_payload_only));
               phase_compensation_lts = exp(1j*angle(curr_all_ltses.*conj_mult_all_ltses));

               rx_payload_and_lts(pkt_idx,global_ant_idx,:) = curr_payload_and_lts.*phase_compensation_payload_lts ;
               rx_payload_only(pkt_idx,global_ant_idx,:) = curr_payload_only.*phase_compensation_payload;
               rx_all_ltses(pkt_idx,global_ant_idx,:) = curr_all_ltses;
               rx_all_ltses_comp(pkt_idx,global_ant_idx,:) = ref_rx_all_ltses(pkt_idx,i,:);
%                else
%                    rx_payload_and_lts(pkt_idx,global_ant_idx,:) = curr_payload_and_lts;
%                    rx_payload_only(pkt_idx,global_ant_idx,:) = curr_payload_only;
%                    rx_all_ltses(pkt_idx,global_ant_idx,:) = curr_all_ltses;
%                end
               
               all_lts = squeeze(rx_all_ltses(pkt_idx,global_ant_idx,:));
               all_lts_comp = squeeze(rx_all_ltses_comp(pkt_idx,global_ant_idx,:));
               
               
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
                        curr_lts = all_lts(lts_start_ind_offseted+half_lts+(lts_idx-1)*num_subc:lts_start_ind_offseted+half_lts+(lts_idx)*num_subc-1);
                        curr_ref_lts = all_lts_comp(lts_start_ind_offseted+half_lts+(lts_idx-1)*num_subc:lts_start_ind_offseted+half_lts+(lts_idx)*num_subc-1);
                        
                        curr_lts_f = fft(curr_lts);
                        curr_ref_lts_f = fft(curr_ref_lts);
                        
                        compensated_lts_f = curr_lts_f.*exp(-1j*(angle(curr_ref_lts_f)));
%                         compensated_lts_f = curr_lts_f.*exp(1j*angle(curr_lts_f.*conj(curr_ref_lts_f)));
%                         compensated_lts_f = curr_lts_f.*exp(1j*angle(conj(curr_ref_lts_f)));
                        
                        debug_lts(pkt_idx,lts_idx,user_idx,global_ant_idx,:) = curr_lts;
                       
                        rx_lts(pkt_idx,lts_idx,user_idx,global_ant_idx,:) = ifft(compensated_lts_f );
                   end

    %                lts_start_ind = lts_start_ind+(ofdm_params.NUM_LTS+0.5)*num_subc+ofdm_params.inter_user_zeros;
        %            xline(lts_start_ind,'m--')
               end
           end
        end
    end
    
    chan_debug = false;
    if(chan_debug)
       lts_user1_across_packets = squeeze(rx_lts(:,1,1,8,:));
       lts_user1_across_packets_nocomp = squeeze(debug_lts(:,1,1,8,:));

       lts_user1_across_packets_f = fft(lts_user1_across_packets,64,2);
       lts_user1_across_packets_nocomp_f = fft(lts_user1_across_packets_nocomp,64,2);

       figure(96)
       subplot(2,1,1)
       plot(1:1:64,abs(fftshift(lts_user1_across_packets_f)))
       hold on
       figure(97)
       subplot(2,1,1)
       plot(1:1:64,abs(fftshift(lts_user1_across_packets_nocomp_f)))
        
       figure(96)
       subplot(2,1,2)
       plot(1:1:64,angle(fftshift(lts_user1_across_packets_f)))
       hold on
       figure(97)
       subplot(2,1,2)
       plot(1:1:64,angle(fftshift(lts_user1_across_packets_nocomp_f)))

       disp("check")
   end
    
end


function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins));
end