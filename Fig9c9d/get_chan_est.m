function [raw_rx_samps,chan_est] = get_chan_est(antenna_state, ofdm_params, samp_fac,...
                                                conj_idx, num_lts_to_use, disp_plots,...
                                                bb_gain, raw_rx_samps)
    
    if(nargin==6)
        disp("Setting gain to 11")
        bb_gain=11;
    elseif(nargin==7)
        raw_rx_samps = get_rx_samps_warp(antenna_state,bb_gain);
    end

    rfc_index = 2;
    conj_index = 4;
    
    num_users = ofdm_params.num_users;
    
    rx_samps = raw_rx_samps;
    rx_samps_rfc = rx_samps(rfc_index,:).';
    for oversamp_ind=1:1:samp_fac
        rx_samples_resamped{oversamp_ind} = delayseq2(rx_samps_rfc(oversamp_ind:samp_fac:end),(oversamp_ind-1)/samp_fac);
    end
    rx_samples_resamped{oversamp_ind+1} = resample(rx_samps(conj_index,:),1,samp_fac);    

    clear op_struct_cell op_struct_conj
    [op_struct_conj, ~] = ofdm_rx(rx_samples_resamped{conj_idx}, ofdm_params);

    if(isempty(op_struct_conj))
        disp("LTS not detected, exiting")
        chan_est = -1;
        return;
    end
    
    lts_inds = op_struct_conj.lts_ind;
    for ind=1:1:samp_fac+1
        [op_struct_cell{ind}, ~] = ofdm_rx(rx_samples_resamped{ind}, ofdm_params, lts_inds);
    end

    op_struct = op_struct_cell;
    start_offset=1;
    [~,~,num_subc,num_packets] = size(op_struct{1}.rx_H_est_vec);
    conj_avg_chans = zeros(num_packets,num_users,samp_fac,num_subc);
%     conj_idx =3;

    for user_idx=1:1:num_users
        if(disp_plots)
            figure(user_idx)
        end
        ref_chan = op_struct{conj_idx}.rx_H_est_vec(user_idx,:,:,start_offset:end);
        ref_mult = exp(-1j*angle(ref_chan));

        for ant_idx=1:1:samp_fac
            curr_chan = op_struct{ant_idx}.rx_H_est_vec(user_idx,:,:,start_offset:end);
            conj_chan = squeeze(curr_chan.*ref_mult);
%                 conj_chan_lts_avg = squeeze(mean(conj_chan,1));
            conj_chan_lts_avg = squeeze(mean(conj_chan(ofdm_params.NUM_LTS-num_lts_to_use+1:ofdm_params.NUM_LTS,:,:),1));
            conj_avg_chans(:,user_idx,ant_idx,:) = conj_chan_lts_avg.';
            if(disp_plots)
                subplot(samp_fac,2,(ant_idx-1)*2+1)
                plot(20*log10(abs(conj_chan_lts_avg)),'Color','b')
                title("Relative channel magnitude: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
                hold on
                ylim([-40,20])
                subplot(samp_fac,2,(ant_idx-1)*2+2)
                plot(unwrap(fftshift(angle(conj_chan_lts_avg(ofdm_params.SC_IND_DATA,:))))*180/pi,'Color','b')
                title("Relative channel phase: Samples "+num2str(ant_idx-1)+" :"+num2str(samp_fac)+":end")
                ylim([-200,200])
                hold on
            end
        end
    end
    chan_est = conj_avg_chans;
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig.', nfft);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins));
end
