%% Parameters and Flags
clear rx_samples rx_samples_resamped op_struct
ofdm_params = wcsng_ofdm_param_gen(64);

% For OFDM RX
count = 1e6; % roughly 40e6 is 1sec data for 40.96MSPS receiver.
test_ids = [0,1,2,3,4,5,6,7,8,9];
% test_ids = [0,1,2,3,4];
test_str="sup_check_after_warp_10m";
start = 0.5e6; % remove initial 1e6 samples because they might be corrupted
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 2;
[tx_samples,ofdm_params] = ofdm_tx(ofdm_params); % Take tx packet and params from OFDM_tx.m
% file_rx_path = "/mnt/ramdisk/";
file_rx_path="/home/wcsng-20/ag_data/"+test_str+"_07_27_2/";
samp_fac = 5;
pow_rf_vec = zeros([numel(test_ids),samp_fac]);
num_cols=2;
%% Loading samples
% Filenames

for test_idx=1:1:numel(test_ids)
    % rx_samples = zeros(count,4);
    % rx_samples_resamped = zeros(count/4),4);

    

    for ind = 1:1
        rx_samples(:,ind) = read_complex_binary(file_rx_path+"rx"+num2str(ind-1)+"_"+test_str+num2str(test_ids(test_idx))+".dat",count,start);
        rx_samples_resamped{ind} = resample(rx_samples(:,ind),1,samp_fac);
    end
    
    for ind=2:2
        rx_samples(:,ind) = read_complex_binary(file_rx_path+"rx"+num2str(ind-1)+"_"+test_str+num2str(test_id)+".dat",count,start);
        for oversamp_ind=1:1:samp_fac
            rx_samples_resamped{1+oversamp_ind} = delayseq2(rx_samples(oversamp_ind:samp_fac:end,ind),(oversamp_ind-1)/samp_fac);
        end
    end

    % figure(1)
    % plot(abs(rx_samples_resamped{1}))

    %% Decode
    ofdm_params.DO_DECODE = 1;
    ofdm_params.DO_APPLY_CFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_SFO_CORRECTION = 1;
    ofdm_params.DO_APPLY_PHASE_ERR_CORRECTION = 1;
    ofdm_params.plot_flag = 0;

    if(ofdm_params.plot_flag==1)
        [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{3}, ofdm_params);
    else 
%         for ind = 1:samp_fac+1
%     %         tic
%             [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params);
%     %         toc
%         end
        for ind = 1:samp_fac+1
            if(ind==1)
                [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params);
            else
                [op_struct(ind), ofdm_params_out(ind)] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params, op_struct(1).lts_ind);
            end
        end
    end
    %%
    start_offset = 5;
%     figure(2)
    
    for ant_idx=2:1:(1+samp_fac)
        conj_idx = 1;
    %     subplot(samp_fac,2,(ant_idx-2)*2+1)
        conj_chan = op_struct(ant_idx).rx_H_est_vec(:,start_offset:end).*conj(op_struct(conj_idx).rx_H_est_vec(:,start_offset:end));
        pow_rf_vec(test_idx,ant_idx-1) = 10*log10(mean(mean(abs(conj_chan(ofdm_params.SC_IND_DATA,:)).^2)));
        %     plot(20*log10(abs(conj_chan)))
    %     title("Relative channel magnitude: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
    %     ylim([-30,10])
    %     subplot(samp_fac,2,(ant_idx-2)*2+2)
    %     plot(unwrap(fftshift(angle(conj_chan(ofdm_params.SC_IND_DATA,:))))*180/pi)
    %     title("Relative channel phase: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
    %     ylim([-230,230])
    end
    figure(1)
    subplot(ceil(numel(test_ids)/num_cols),num_cols,num_cols*floor(test_idx/num_cols)+mod(test_idx,num_cols))
    plot(pow_rf_vec(test_idx,:),'Marker','*','Markersize',20)
%     ylim([-10,0])
end
%%
% figure(3)
% for i=1:1:1
%     subplot(1,2,i)
%     plot(20*log10(abs(op_struct(i).rx_H_est_vec(:,start_offset:end))))
% end
function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins));
end