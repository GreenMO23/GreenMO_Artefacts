%% Parameters and Flags
clear op_struct rx_samps raw_rx_samps rx_samples_resamped
ofdm_params = wcsng_ofdm_param_gen(64);

% For OFDM RX
count = 2e6; % roughly 40e6 is 1sec data for 40.96MSPS receiver.
test_id = 3 ;
start = 0.5e6; % remove initial 1e6 samples because they might be corrupted
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 2;
[tx_samples,ofdm_params] = ofdm_tx(ofdm_params); % Take tx packet and params from OFDM_tx.m

mat_source = true;

samp_fac = 4;
start_offset = 1;
% warp_expmt_base_str = "rx_samps_warp_40Mswitch_hg_14";

savpath_dir = "/home/wcsng-20/ag_data/MMSE_baseline_expmt";
expmt_str="all_four_ants_rot";
path_dir = savpath_dir+"/"+expmt_str+"/";
expmt_mode = "pcb";

% test_ids = [0,1,2,3,4]+5;
conj_index = 1;
rfc_index = 2;
num_tests = 6;
test_ids = 0:1:num_tests-1;
test_off = 0;
test_ids = test_ids+test_off;
% test_ids = [0,1,2,3,4];
num_cols=2;
for test_idx=1:1:numel(test_ids)
    fname = path_dir+expmt_mode+"_"+num2str(test_ids(test_idx));
    load(fname);
    %% Downsample and phase align appropriately
    if(expmt_mode=="pcb")
        rx_samples_resamped{1} = resample(rx_samps(conj_index,:),1,samp_fac);    
        rx_samps_rfc = rx_samps(rfc_index,:).';
        for oversamp_ind=1:1:samp_fac
            rx_samples_resamped{1+oversamp_ind} = delayseq2(rx_samps_rfc(oversamp_ind:samp_fac:end),(oversamp_ind-1)/samp_fac);
        end

    else
        for ind=1:samp_fac  
            rx_samples_resamped{ind} = resample(rx_samps(ind,:),1,samp_fac);
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
    figure(3)
    subplot(ceil(numel(test_ids)/num_cols),num_cols,num_cols*floor(test_idx/num_cols)+mod(test_idx,num_cols))
    pow_rf_vec = zeros([1,samp_fac]);
    if(expmt_mode=="pcb")
        num_rfc_resamp = 1+samp_fac;
        start_ind = 2;
    else
        num_rfc_resamp = samp_fac;
        start_ind = 1;
    end
    for ant_idx=start_ind:1:num_rfc_resamp
            
        %     subplot(samp_fac,2,(ant_idx-2)*2+1)
            conj_chan = op_struct(ant_idx).rx_H_est_vec(:,start_offset:end).*conj(op_struct(conj_index).rx_H_est_vec(:,start_offset:end));
            pow_rf_vec(ant_idx-(start_ind-1)) = 10*log10(mean(mean(abs(conj_chan(ofdm_params.SC_IND_DATA,:)).^2)));
            %     plot(20*log10(abs(conj_chan)))
        %     title("Relative channel magnitude: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
        %     ylim([-30,10])
        %     subplot(samp_fac,2,(ant_idx-2)*2+2)
        %     plot(unwrap(fftshift(angle(conj_chan(ofdm_params.SC_IND_DATA,:))))*180/pi)
        %     title("Relative channel phase: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
        %     ylim([-230,230])
    end

    plot(pow_rf_vec,'Marker','*','Markersize',20)
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins));
end
