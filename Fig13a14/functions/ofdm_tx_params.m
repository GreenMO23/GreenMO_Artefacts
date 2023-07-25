function [tx_waveforms,loopback,ofdm_tx_params_allusers] = ofdm_tx_params(meta_params,hwr_flag)
    clear ofdm_tx_params_allusers tx_waveforms loopback

    plot_flag = meta_params.plot_setting;
    
    if(nargin==1)
        hwr_flag=false;
    end
    
    num_users = meta_params.num_users;
    
    ofdm_tx_params = wcsng_ofdm_param_gen(64);
    if(hwr_flag)
        num_packets = 400;
    else
        num_packets = 2;
    end
    
    ofdm_tx_params.packet_step_size = 1;
    ofdm_tx_params.MOD_ORDER = 16;
    ofdm_tx_params.channel_coding_rate = 0.75;
    ofdm_tx_params.NUM_LTS = 10;
    ofdm_tx_params.num_users = num_users;
    ofdm_tx_params.enable_channel_codes=true;
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

    if(plot_flag)
        figure(2)
    end
    
    

    for i=1:1:ofdm_tx_params.num_users
        ofdm_tx_params.user_index = i-1; % required by current ofdm_tx class to separate the LTS's
        if(hwr_flag)
            filename_tx = path_dir+"uhd_ofdm_tx_"+num2str(ofdm_tx_params.num_users)+"user_coded_QAM16_34r_"+num2str(i-1)+".dat";
            [tx_waveforms(i),loopback(i),ofdm_tx_params_allusers(i)] = ofdm_tx(ofdm_tx_params); % Take tx packet from OFDM_tx.m
            
            
            %     rms(tx_samples)

            %     tx_packet_tosave = repmat(tx_samples,num_packets,1);
            %     size(tx_packet_tosave)
            %     write_complex_binary(tx_packet_tosave,filename_tx);
            %     clear tx_packet_tosave
        else
            [tx_sig,loopback_sig, ofdm_tx_params_allusers{i}]= ofdm_tx(ofdm_tx_params); % Take tx packet from OFDM_tx.m;
            sig_pow = 10*log10(rms(tx_sig).^2/0.001); % dbm
            mult_fac = sqrt(10^((meta_params.tx_pow-sig_pow)/10));
            tx_samples = mult_fac*tx_sig;
            loopback_samples = mult_fac*loopback_sig;
            tx_waveforms{i} = repmat(tx_samples,num_packets,1);
            loopback{i} = repmat(loopback_samples,num_packets,1);
            num_zeros_append = 128;
            tx_waveforms{i} = [zeros(1,num_zeros_append) tx_waveforms{i}.' zeros(1,num_zeros_append)].';
            loopback{i} = [zeros(1,num_zeros_append) loopback{i}.' zeros(1,num_zeros_append)].';
        end
        
        if(plot_flag)
            plot(abs(tx_waveforms{i}))
            hold on
        end
    end

end