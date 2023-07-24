%% Parameters and Flags
clear op_struct rx_samps raw_rx_samps rx_samples_resamped ofdm_params_out
ofdm_params = wcsng_ofdm_param_gen(64);

% For OFDM RX
count = 2e6; % roughly 40e6 is 1sec data for 40.96MSPS receiver.
test_id = 1 ;
start = 0.5e6; % remove initial 1e6 samples because they might be corrupted
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 4;
ofdm_params.NUM_LTS = 10;
num_lts_to_use = 2;
ofdm_params.num_users = 4;
num_users = ofdm_params.num_users;
[tx_samples,ofdm_params] = ofdm_tx(ofdm_params); % Take tx packet and params from OFDM_tx.m

mat_source = true;
nrf = 2;
samp_fac = 4;
savpath_dir = "/home/wcsng-20/ag_data/MMSE_baseline_expmt_decode_4user";
expmt_str="0";
path_dir = savpath_dir+"/"+expmt_str+"/";
expmt_mode = "mmse";
% fname = savpath_dir+expmt_str+"_pcb_"+num2str(test_ids(test_idx));
test_idx=0;
fname = path_dir+expmt_mode+"_"+num2str(test_idx);
load(fname);

if(expmt_mode=="pcb")
    sum_idx=1;
    h1_idx= 2;
    h2_idx= 3;
    zero_idx = 4;
    sup_plot=false;
else
    sup_plot=false;
end


% expmt_mode = "pcb";

%% Downsample and phase align appropriately
if(expmt_mode=="pcb")
        rx_samps_rfc = rx_samps(rfc_index,:).';
        for oversamp_ind=1:1:samp_fac
            rx_samples_resamped{oversamp_ind} = delayseq2(rx_samps_rfc(oversamp_ind:samp_fac:end),(oversamp_ind-1)/samp_fac);
        end
        rx_samples_resamped{oversamp_ind+1} = resample(rx_samps(conj_index,:),1,samp_fac);    

else
    for ind=1:samp_fac  
        rx_samples_resamped{ind} = resample(rx_samps(ind,:),1,samp_fac);
%            rx_samples_resamped{ind} = rx_samps(ind,1:samp_fac:end);
    end
end
%% decode
ofdm_params.DO_DECODE = 1;
ofdm_params.DO_APPLY_CFO_CORRECTION = 1;
ofdm_params.DO_APPLY_SFO_CORRECTION = 1;
ofdm_params.DO_APPLY_PHASE_ERR_CORRECTION = 1;
ofdm_params.plot_flag = 0;

if(ofdm_params.plot_flag==1)
    [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{3}, ofdm_params);
else 
    if(expmt_mode=="pcb")
        num_rfc_resamp = 1+samp_fac;
    else
        num_rfc_resamp = samp_fac;
    end

    for ind = 1:num_rfc_resamp
        tic
        if(ind==1)
            [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params);
            lts_inds = op_struct(1).lts_ind;
        else
            [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params, op_struct(1).lts_ind);
        end
        toc
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
    figure(user_idx)
    ref_chan = op_struct(conj_idx).rx_H_est_vec(user_idx,:,:,start_offset:end);
    ref_mult = exp(-1j*angle(ref_chan));
    for ant_idx=1:1:samp_fac
        curr_chan = op_struct(ant_idx).rx_H_est_vec(user_idx,:,:,start_offset:end);
        conj_chan = squeeze(curr_chan.*ref_mult);
        conj_chan_lts_avg = squeeze(mean(conj_chan(ofdm_params.NUM_LTS-num_lts_to_use+1:ofdm_params.NUM_LTS,:,:),1));
        conj_avg_chans(:,user_idx,ant_idx,:) = conj_chan_lts_avg.';
        
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
    
%%
figure(3)
for i=1:1:4
    subplot(2,2,i)
    plot(abs(rx_samples_resamped{i}))
    hold on
end
lts_len = num_subc;
payload_len = (ofdm_params.N_SC+ofdm_params.CP_LEN)*ofdm_params.N_OFDM_SYMS;
rx_lts = zeros(num_packets,ofdm_params.NUM_LTS,num_users,samp_fac,num_subc);
rx_payload = zeros(num_packets, samp_fac,payload_len);
for pkt_idx=1:1:num_packets
   for i=1:1:4
       pkt_lts_ind = op_struct(i).packet_lts_list(pkt_idx);
       pkt_payload_ind = op_struct(i).packet_payload_list(pkt_idx);
       raw_samps_rfi = rx_samples_resamped{i};
       rx_payload(pkt_idx,i,:) = raw_samps_rfi(pkt_payload_ind:pkt_payload_ind+payload_len-1);
       subplot(2,2,i)
       xline(pkt_lts_ind,'r--')
       xline(pkt_payload_ind,'g-.')
       xline(pkt_payload_ind+payload_len,'m:')
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
nonnull_subc = [2:27 39:64];
num_nonnull_subc = numel(nonnull_subc);
combined_lts_f = zeros(num_packets,num_users,num_users,ofdm_params.NUM_LTS,num_subc);
user_comb_vecs = zeros(num_users,samp_fac,num_subc);
for pkt_idx=1:1:num_packets
    rx_lts_f = fft(squeeze(rx_lts(pkt_idx,:,:,:,:)),num_subc,4);% 4th axis is subcarrier
    for subc_idx=1:1:num_nonnull_subc 
        chan_mat = chan_to_use(:,:,nonnull_subc(subc_idx));
        user_vec = 1:1:num_users;
        for user_idx=1:1:num_users
            curr_user = user_idx;
            interf_users = user_vec;
            interf_users(curr_user) =[];
            null_chan = null(chan_mat(interf_users,:));
            sig_chan = chan_mat(curr_user,:);
            
            null_comb_vec = sig_chan*null_chan;
            fin_comb_vec = null_chan*null_comb_vec';
            
            user_comb_vecs(user_idx,:,nonnull_subc(subc_idx)) = fin_comb_vec;
            
            for comb_user_idx=1:1:num_users
                curr_lts_to_be_combined = squeeze(rx_lts_f(:,comb_user_idx,:,nonnull_subc(subc_idx)));
                repmat_combvec = repmat(fin_comb_vec.',[ofdm_params.NUM_LTS,1]);
                combined_samp = sum(curr_lts_to_be_combined.*repmat_combvec,2);
                combined_lts_f(pkt_idx,user_idx,comb_user_idx,:,nonnull_subc(subc_idx))=combined_samp;     
            end
        end
    end
end

% pkt,user,user,10lts,nonnulls
combined_lts_f_non_null = combined_lts_f(:,:,:,:,nonnull_subc);
sinr_calc = zeros(num_packets,num_users,num_nonnull_subc );
for pkt_idx=1:1:num_packets
    for user_idx=1:1:num_users
        subc_pow_profile = (abs(squeeze(combined_lts_f_non_null(pkt_idx,user_idx,:,:,:))).^2);
        mean_pow_per_comb_user = squeeze(mean(subc_pow_profile,2)); %mean across 10lts
        noise_pow = sqrt(squeeze(var(subc_pow_profile(user_idx,:,:),1)));
        sig_pow = mean_pow_per_comb_user(user_idx,:).';
        curr_user = user_idx;
        interf_users = user_vec;
        interf_users(curr_user) =[];
        interf_pow = sum(mean_pow_per_comb_user(interf_users,:),1).';
        sinr_vec = sig_pow./(interf_pow+noise_pow);
        sinr_calc(pkt_idx,user_idx,:) = sinr_vec;
        
    end
end

mean(mean(mean(sinr_calc)))
%% Across SINR
figure(4)
pkt_idx=5;
scatter(1:1:num_nonnull_subc, squeeze(sinr_calc(pkt_idx,1,:)))
hold on
scatter(1:1:num_nonnull_subc, squeeze(sinr_calc(pkt_idx,2,:)))
legend("User 1 "+expmt_mode, "User 2 "+expmt_mode )

%% Across packets
figure(5)
sinr_avg_subc = squeeze(mean(sinr_calc,3));
plot(1:1:num_packets,sinr_avg_subc(:,1))
hold on
plot(1:1:num_packets,sinr_avg_subc(:,2))
legend("User 1 "+expmt_mode, "User 2 "+expmt_mode )
ylim([0,50])
% mean(mean(mean(sinr_calc)))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug Codes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pow_rfc_debug=false;
if(pow_rfc_debug)
    for user_idx=1:1:num_users
        figure(2+user_idx)
        ref_chan = op_struct(conj_idx).rx_H_est_vec(user_idx,:,:,start_offset:end);
        ref_mult = exp(-1j*angle(ref_chan));
        pow_rf_vec = zeros([1,samp_fac]);
        for ant_idx=start_ind:1:num_rfc_resamp

                curr_chan = op_struct(ant_idx).rx_H_est_vec(user_idx,:,:,start_offset:end);
                conj_chan = squeeze(curr_chan.*ref_mult);
                conj_chan_lts_avg = squeeze(mean(conj_chan,1));

                pow_rf_vec(ant_idx-(start_ind-1)) = 10*log10(mean(mean(abs(conj_chan_lts_avg(ofdm_params.SC_IND_DATA,:)).^2)));

        end
        plot(pow_rf_vec,'Marker','*','Markersize',20)
    end
end



%%
ant_debug=false;
if(ant_debug)
    for ant_idx=1:1:1+samp_fac
        subplot(samp_fac+1,2,(ant_idx-1)*2+1)
        conj_chan = op_struct(ant_idx).rx_H_est_vec(:,start_offset:end);
        plot(20*log10(abs(conj_chan)))

        title("Relative channel magnitude: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
        ylim([-20,30])
        subplot(samp_fac+1,2,(ant_idx-1)*2+2)
        plot(unwrap(fftshift(angle(conj_chan(ofdm_params.SC_IND_DATA,:))))*180/pi)
        title("Relative channel phase: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")

        ylim([-300,300])
    end
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins));
end

