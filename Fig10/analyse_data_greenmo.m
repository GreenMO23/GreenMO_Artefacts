function [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
           avg_sinr_users_across_packets,min_sinr_users_across_packets,...
           avg_ber_users_across_packets,max_ber_users_across_packets]...
           = analyse_data_greenmo(posn_index,num_ant,data_dir, test_ids)

clear op_struct op_struct_cell rx_samps raw_rx_samps rx_samples_resamped ofdm_params_out 
clear avg_sinr avg_ber avg_snr_evm ber_dict evm_snr_dict sinr_dict
                    
%% Parameters and Flags
expmt_mode = "bit_phase";                    
ofdm_params = wcsng_ofdm_param_gen(64);
ofdm_params.num_users = 4;

% For OFDM RX
count = 2e6; % roughly 40e6 is 1sec data for 40.96MSPS receiver.
test_id = 1 ;
start = 0.5e6; % remove initial 1e6 samples because they might be corrupted
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 16;
ofdm_params.NUM_LTS = 10;
num_lts_to_use = 10;
ofdm_params.inter_user_zeros = 512;

num_users = ofdm_params.num_users;
ofdm_params.enable_channel_codes = true;



savpath_dir = data_dir;
tx_params_dir = savpath_dir+"/tx_data"+num2str(num_users)+"/";
load(tx_params_dir+"ofdm_params_tx.mat");



expmt_str="rx_data"+num2str(num_ant)+"_"+num2str(posn_index);

path_dir = savpath_dir+"/"+expmt_str+"/";



[tx_samples,ofdm_params] = ofdm_tx(ofdm_params); % Take tx packet and params from OFDM_tx.m
plot_debug=false;
disp_plots=false;
nrf = 2;
samp_fac = 4;


num_tests = numel(test_ids);
avg_evm_snr_users_across_packets = zeros([num_tests,num_users]);
min_evm_snr_users_across_packets = zeros([num_tests,num_users]);
avg_sinr_users_across_packets = zeros([num_tests,num_users]);
min_sinr_users_across_packets = zeros([num_tests,num_users]);
avg_ber_users_across_packets = zeros([num_tests,num_users]);
max_ber_users_across_packets = zeros([num_tests,num_users]);


for test_index = 1:1:num_tests
    % Load data
    test_idx = test_ids(test_index);
%     disp("Expmt: "+expmt_mode+", Test idx: "+num2str(test_idx)+", Users: "+num2str(num_users));
    fname = path_dir+"rx_samps_"+expmt_mode+"_"+num2str(test_idx);
    load(fname);
    rx_samps = raw_rx_samps;


    % expmt_mode = "pcb";

    %% Downsample and phase align appropriately
    rfc_index = 2;
    conj_index = 4;
    rx_samps_rfc = rx_samps(rfc_index,:).';
    for oversamp_ind=1:1:samp_fac
        rx_samples_resamped{oversamp_ind} = delayseq2(rx_samps_rfc(oversamp_ind:samp_fac:end),(oversamp_ind-1)/samp_fac);
    end
    rx_samples_resamped{oversamp_ind+1} = resample(rx_samps(conj_index,:),1,samp_fac);    

    %% decode
    ofdm_params.DO_DECODE = 1;
    ofdm_params.DO_APPLY_CFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_SFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_PHASE_ERR_CORRECTION = 1;

    break_for = false;
    for ind = 1:1:samp_fac
%                 tic
        if(ind==1)
            [op_struct_cell{ind}, ofdm_params_out{ind}] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params);
            if(isempty(op_struct_cell{ind}))
                %  % what we need for each expmt  
                avg_sinr(test_index) = -1;
                avg_snr_evm(test_index) = -1;
                avg_ber(test_index) = -1;

                ber_dict{test_index} = -1;
                evm_snr_dict{test_index} = -1;
                sinr_dict{test_index} = -1;
                break_for = true;
                break;
            else
                lts_inds = op_struct_cell{1}.lts_ind;
            end
        else
            [op_struct_cell{ind}, ofdm_params_out{ind}] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params, lts_inds);
        end
%                 toc

    end

    if(break_for)
        continue
    else
        for ind=1:1:samp_fac
            op_struct(ind) =op_struct_cell{ind};
        end
    end
    %%
    start_offset = 1;
    % if(expmt_mode=="pcb")
    %     num_rfc_resamp = 1+samp_fac;
    % else
    %     num_rfc_resamp = samp_fac;
    % end
    conj_idx = 1;
    start_offset = 1;
    [~,~,num_subc,num_packets] = size(op_struct(conj_idx).rx_H_est_vec);
    conj_avg_chans = zeros(num_packets,num_users,samp_fac,num_subc);

    for user_idx=1:1:num_users
        if(plot_debug)
            figure(user_idx)
        end
        ref_chan = op_struct(conj_idx).rx_H_est_vec(user_idx,:,:,start_offset:end);
        ref_mult = exp(-1j*angle(ref_chan));
        for ant_idx=1:1:samp_fac
            curr_chan = op_struct(ant_idx).rx_H_est_vec(user_idx,:,:,start_offset:end);
            conj_chan = squeeze(curr_chan.*ref_mult);
%                 conj_chan_lts_avg = squeeze(mean(conj_chan,1));
            conj_chan_lts_avg = squeeze(mean(conj_chan(ofdm_params.NUM_LTS-num_lts_to_use+1:ofdm_params.NUM_LTS,:,:),1));
            conj_avg_chans(:,user_idx,ant_idx,:) = conj_chan_lts_avg.';
            if(plot_debug)
                subplot(samp_fac,2,(ant_idx-1)*2+1)
                plot(20*log10(abs(conj_chan_lts_avg)))

                if(ant_idx==sum_idx+1 && sup_plot)
                    hold on
                    conj_chan_A = op_struct(h1_idx).rx_H_est_vec(user_idx,:,:,start_offset:end).*ref_mult;
                    conj_chan_B = op_struct(h2_idx).rx_H_est_vec(user_idx,:,:,start_offset:end).*ref_mult;
                    zero_chan =  op_struct(zero_idx).rx_H_est_vec(user_idx,:,:,start_offset:end).*ref_mult;
                    sum_chan = (mean(conj_chan_A,1)+mean(conj_chan_B,1))-mean(zero_chan,1);
                    plot(20*log10(abs(sum_chan)))
                end

                title("Relative channel magnitude: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
                ylim([-20,35])
                subplot(samp_fac,2,(ant_idx-1)*2+2)
                plot(unwrap(fftshift(angle(conj_chan_lts_avg(ofdm_params.SC_IND_DATA,:))))*180/pi)
                title("Relative channel phase: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")

                if(ant_idx==sum_idx+1 && sup_plot)
                    hold on
                    plot(unwrap(fftshift(angle(sum_chan(ofdm_params.SC_IND_DATA,:))))*180/pi);
                end
                ylim([-300,300])
            end
        end
    end

    %%
%         plot_debug=true
    if(plot_debug)
        figure(3)
        for i=1:1:4
            subplot(2,2,i)
            plot(abs(rx_samples_resamped{i}))
            hold on
        end
    end
    sing_lts_len = num_subc;
    lts_len = (((ofdm_params.NUM_LTS+0.5)*sing_lts_len)+ofdm_params.inter_user_zeros)*num_users;
    payload_len = (ofdm_params.N_SC+ofdm_params.CP_LEN)*ofdm_params.N_OFDM_SYMS;
    pkt_and_payload_plus_zeros = (ceil((lts_len+payload_len)/ofdm_params.N_SC)+1)*ofdm_params.N_SC;
    num_subc_groups = pkt_and_payload_plus_zeros/ofdm_params.N_SC;
    rx_lts = zeros(num_packets,ofdm_params.NUM_LTS,num_users,samp_fac,num_subc);
    rx_payload_and_lts = zeros(num_packets, samp_fac,pkt_and_payload_plus_zeros);
    for pkt_idx=1:1:num_packets
       for i=1:1:4
           pkt_lts_ind = op_struct(i).packet_lts_list(pkt_idx);
           pkt_payload_ind = op_struct(i).packet_payload_list(pkt_idx);
           raw_samps_rfi = rx_samples_resamped{i};
           rx_payload_and_lts(pkt_idx,i,:) = raw_samps_rfi(pkt_lts_ind:pkt_lts_ind+pkt_and_payload_plus_zeros-1);
           if(plot_debug)
               subplot(2,2,i)
               xline(pkt_lts_ind,'r--')
               xline(pkt_payload_ind,'g-.')
               xline(pkt_lts_ind+pkt_and_payload_plus_zeros,'m:')
           end
           all_lts = raw_samps_rfi(pkt_lts_ind:pkt_payload_ind-1);
           lts_offset=0;
           half_lts = num_subc/2;
    %        figure(10)
    %        plot(abs(all_lts))
    %        hold on
           lts_start_ind = 1;
           for user_idx=1:1:num_users
               for lts_idx=1:1:ofdm_params.NUM_LTS
    %                xline(lts_start_ind+half_lts+(lts_idx-1)*num_subc,'r:') 
    %                 (num_packets,ofdm_params.NUM_LTS,num_users,samp_fac,num_subc)
                   rx_lts(pkt_idx,lts_idx,user_idx,i,:) = all_lts(lts_start_ind+half_lts+(lts_idx-1)*num_subc:lts_start_ind+half_lts+(lts_idx)*num_subc-1);
               end
    %            xline(lts_start_ind+half_lts+(lts_idx)*num_subc,'r:') 
               lts_start_ind = lts_start_ind+(ofdm_params.NUM_LTS+0.5)*num_subc+ofdm_params.inter_user_zeros;
    %            xline(lts_start_ind,'m--')
           end
       end
    end
    %% SINR calculations now finally
    % conj_avg_chans = zeros(num_packets,num_users,samp_fac,num_subc);
    % rx_lts = zeros(num_packets,ofdm_params.NUM_LTS,num_users,samp_fac,num_subc);
    pkt_to_use = 1;
    chan_to_use = squeeze(conj_avg_chans(pkt_to_use,:,:,:));
%         chan_to_use = squeeze(mean(conj_avg_chans,[1]));

    nonnull_subc = [2:27 39:64];
    num_nonnull_subc = numel(nonnull_subc);
    combined_lts_f = zeros(num_packets,num_users,num_users,ofdm_params.NUM_LTS,num_subc);
    user_comb_vecs = zeros(num_users,samp_fac,num_subc);
    tdma_snrs = zeros(num_packets,num_users,samp_fac,num_nonnull_subc);
    % Take FFTCombine and take IFFT to get back to td
    % rx_payload_and_lts = zeros(num_packets, samp_fac,pkt_and_payload_plus_zeros);
    rx_payload_and_lts_groups_64 = reshape(rx_payload_and_lts,num_packets,samp_fac,num_subc,num_subc_groups);
    rx_payload_and_lts_groups_64_f = fft(rx_payload_and_lts_groups_64,64,3); 
    combined_payload_and_lts_f = zeros(num_packets,num_users,num_subc,num_subc_groups);


    for pkt_idx=1:1:num_packets
        if(num_users>1)
            rx_lts_f = fft(squeeze(rx_lts(pkt_idx,:,:,:,:)),num_subc,4);% 4th axis is subcarrier
        else
            rx_lts_f = fft(squeeze(rx_lts(pkt_idx,:,:,:,:)),num_subc,3);% 4th axis is subcarrier
        end
        for subc_idx=1:1:num_nonnull_subc 
            if(num_users==1)
                chan_mat = chan_to_use(:,nonnull_subc(subc_idx));
            else
                chan_mat = chan_to_use(:,:,nonnull_subc(subc_idx));
            end
            user_vec = 1:1:num_users;
            for user_idx=1:1:num_users
                if(num_users>1)
                    curr_user = user_idx;
                    interf_users = user_vec;
                    interf_users(curr_user) =[];
                    null_chan = null(chan_mat(interf_users,:));
                    sig_chan = chan_mat(curr_user,:);

                    null_comb_vec = sig_chan*null_chan;
                    fin_comb_vec = null_chan*null_comb_vec';

                    user_comb_vecs(user_idx,:,nonnull_subc(subc_idx)) = fin_comb_vec;

                    curr_payload_to_be_combined = squeeze(rx_payload_and_lts_groups_64_f(pkt_idx,:,nonnull_subc(subc_idx),:)).';
                    repmat_combvec = repmat(fin_comb_vec.',[num_subc_groups,1]);
                    combined_samp = sum(curr_payload_to_be_combined.*repmat_combvec,2);
                    combined_payload_and_lts_f(pkt_idx,user_idx,nonnull_subc(subc_idx),:)=combined_samp;
                    curr_lts_tdmasnr = squeeze(rx_lts_f(:,user_idx,:,nonnull_subc(subc_idx)));

%                         sig_pow = squeeze(mean(abs(curr_lts_to_be_combined).^2,1));
%                         sig_pow = squeeze(mean(abs(curr_lts_to_be_combined),1));
%                         noise_pow = sqrt(squeeze(var(abs(curr_lts_to_be_combined).^2,1)));
                    sig_pow = squeeze(mean(abs(curr_lts_tdmasnr).^2,1));
                    noise_pow = squeeze(var(abs(curr_lts_tdmasnr),1));
                    tdma_snrs(pkt_idx,user_idx,:,subc_idx) = sig_pow./noise_pow;

                    for comb_user_idx=1:1:num_users
                        curr_lts_to_be_combined = squeeze(rx_lts_f(:,comb_user_idx,:,nonnull_subc(subc_idx)));

%                         sig_pow = squeeze(mean(abs(curr_lts_to_be_combined).^2,1));
%                         sig_pow = squeeze(mean(abs(curr_lts_to_be_combined),1));
%                         noise_pow = sqrt(squeeze(var(abs(curr_lts_to_be_combined).^2,1)));
%                             sig_pow = squeeze(mean(abs(curr_lts_to_be_combined).^2,1));
%                             noise_pow = squeeze(var(abs(curr_lts_to_be_combined),1));
%                             tdma_snrs(pkt_idx,user_idx,:,subc_idx) = sig_pow./noise_pow;

                        repmat_combvec = repmat(fin_comb_vec.',[ofdm_params.NUM_LTS,1]);
                        combined_samp = sum(curr_lts_to_be_combined.*repmat_combvec,2);
                        combined_lts_f(pkt_idx,user_idx,comb_user_idx,:,nonnull_subc(subc_idx))=combined_samp;     
                    end
                else
                    ref_chain_idx =1;
                    curr_payload_to_be_combined = squeeze(rx_payload_and_lts_groups_64_f(pkt_idx,ref_chain_idx,nonnull_subc(subc_idx),:)).';
                    curr_lts_to_be_combined = squeeze(rx_lts_f(:,ref_chain_idx,nonnull_subc(subc_idx)));

                    sig_pow = squeeze(mean(abs(curr_lts_to_be_combined).^2,1));
                    noise_pow = squeeze(var(abs(curr_lts_to_be_combined),1));
                    tdma_snrs(pkt_idx,user_idx,:,subc_idx) = sig_pow./noise_pow;

                    combined_lts_f(pkt_idx,user_idx,user_idx,:,nonnull_subc(subc_idx))=curr_lts_to_be_combined;
                    combined_payload_and_lts_f(pkt_idx,user_idx,nonnull_subc(subc_idx),:)=curr_payload_to_be_combined;
                end
            end
        end
    end

    combined_payload_and_lts_t = ifft(combined_payload_and_lts_f,64,3);
    combined_payload_and_lts = reshape(combined_payload_and_lts_t,num_packets,num_users,pkt_and_payload_plus_zeros,[]);
    if(plot_debug)
        pkt_idx = 1;
        figure(8)
        num_rows=2;
        num_cols=ceil(num_users/num_rows); 
        for i=1:1:num_users
            subplot(num_cols,num_rows,i)
            plot(squeeze(abs(combined_payload_and_lts(1,i,:))))
%                 hold on
        end
    end
    % pkt,user,user,10lts,nonnulls
    combined_lts_f_non_null = combined_lts_f(:,:,:,:,nonnull_subc);
    sinr_calc = zeros(num_packets,num_users,num_nonnull_subc );
    ber_calc = zeros(num_packets,num_users);
%         evm_calc = zeros(num_packets,num_users);
    evm_snr_calc = zeros(num_packets,num_users);
    clear all_users_decode_struct
    for pkt_idx=1:1:num_packets
        clear decode_params decode_op_structs
        for user_idx=1:1:num_users
            subc_pow_profile = (abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:))).^2);
            if(num_users>1)
                mean_pow_per_comb_user = squeeze(mean(subc_pow_profile,2)); %mean across 10lts
                subc_mag_profile = abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:)));
                noise_pow = squeeze(var(subc_mag_profile(user_idx,:,:),1));
                sig_pow = mean_pow_per_comb_user(user_idx,:).';
                curr_user = user_idx;
                interf_users = user_vec;
                interf_users(curr_user) =[];
                interf_pow = sum(mean_pow_per_comb_user(interf_users,:),1).';
                sinr_vec = sig_pow./(interf_pow+noise_pow);
                sinr_calc(pkt_idx,user_idx,:) = sinr_vec;
            else
                mean_pow_per_comb_user = squeeze(mean(subc_pow_profile,1)); %mean across 10lts
                subc_mag_profile = abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:)));
                noise_pow = squeeze(var(subc_mag_profile,1));
                sig_pow = mean_pow_per_comb_user(user_idx,:);
                sinr_vec = sig_pow./(noise_pow);
                sinr_calc(pkt_idx,user_idx,:) = sinr_vec;
            end    


            curr_payload = combined_payload_and_lts(pkt_idx,user_idx,:);

            [decode_op_structs(user_idx),decode_params(user_idx)]= ofdm_decode(squeeze(curr_payload).', ofdm_params_tx(user_idx), ofdm_params,user_idx-1);
            ber_calc(pkt_idx, user_idx) = decode_op_structs(user_idx).berratio;
            evm_snr_calc(pkt_idx, user_idx) = decode_op_structs(user_idx).snr;
        end
        all_users_decode_struct{pkt_idx} = decode_op_structs;
    end

%       % what we need for each expmt  
%     avg_sinr_users(test_index,:) = 10*log10(mean(sinr_calc,[1,3])); 
%     avg_sinr(test_index) = 10*log10(mean(mean(mean(sinr_calc))));
%     avg_snr_evm(test_index) = 10*log10(mean(mean(evm_snr_calc)));
%     avg_ber(test_index) = mean(mean(mean(ber_calc)));

%       % what we need for each expmt
    ber_dict{test_index} = ber_calc;
    evm_snr_dict{test_index} = evm_snr_calc;
    sinr_dict{test_index} = sinr_calc;

    avg_evm_snr_users_across_packets(test_index,:) = 10*log10(mean(evm_snr_calc,1));
    min_evm_snr_users_across_packets(test_index,:) = 10*log10(min(evm_snr_calc,[],1));
    
    sinr_subc_averaged = squeeze(mean(sinr_calc,[3])); 
    avg_sinr_users_across_packets(test_index,:) = 10*log10(mean(sinr_subc_averaged,1));
    min_sinr_users_across_packets(test_index,:) = 10*log10(min(sinr_subc_averaged,[],1));
    
    avg_ber_users_across_packets(test_index,:) = mean(ber_calc,1);
    max_ber_users_across_packets(test_index,:) = min(ber_calc,[],1);
    

end
% 
% avg_ber
% avg_sinr
% avg_snr_evm
% path_dir


% avg_tdma_snr = 10*log10(mean(tdma_snrs,[1,3,4]))
% max_user_sinr = squeeze(max(avg_sinr_users,[],2)).'
% capacity_tdma = sum(avg_tdma_snr)/num_users; % defined as log10 for now, correct is log2
% capacity_mimo = sum(max_user_sinr); % defined as log10 for now, correct is log2

% gain_fac = capacity_mimo/capacity_tdma

% if(num_ant==4)

save(path_dir+"/decode_out_fig10.mat",'ber_dict','evm_snr_dict','sinr_dict','tdma_snrs');




end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins));
end

