root_dir = "/media/wcsng-20/57810C9C5890F0C7/";
savpath_dir = root_dir+"ag_data2/8ant_sup_working";
num_users = 1;
dir_prefix = "rx_data99";
custom_suffix = true;

if(~custom_suffix)
    posn_index = 0;
    expmt_str="rx_data"+num2str(num_users)+"_"+num2str(posn_index);
    ret_val = system("mkdir "+savpath_dir+"/"+expmt_str);
    path_dir = savpath_dir+"/"+expmt_str+"/";
    expmt_modes = ["chan_est", "analog_choice","random"];
%     expmt_ids=[0,1,2];
else
    expmt_str="sup_check_0";
%     ret_val = system("mkdir "+savpath_dir+"/"+expmt_str);
    path_dir = savpath_dir+"/"+expmt_str+"/";
    if ~exist(path_dir, 'dir')
       mkdir(path_dir)
    end
    
    rx_samp_file_suffix = "first_data_1";
    
    expmt_modes = ["chan_est"];
    expmt_ids=[0];
end

%% set up ofdm params
ofdm_params = wcsng_ofdm_param_gen(64);
ofdm_params.num_users = 1;
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 16;
ofdm_params.NUM_LTS = 10;
num_lts_to_use = 10;
ofdm_params.inter_user_zeros = 512;
num_users = ofdm_params.num_users;
ofdm_params.enable_channel_codes = true;
tx_params_dir = savpath_dir+"/tx_data"+num2str(num_users)+"/";
load(tx_params_dir+"ofdm_params_tx.mat");
[~,ofdm_params] = ofdm_tx(ofdm_params); % Take tx packet and params from OFDM_tx.m
disp_plots=true;
samp_fac = 4;
user_idx = 1;
%% final experiment loops
clear ch_est_rx_samps 
for expmt_mode=expmt_modes
    for expmt_id=1:1:numel(expmt_ids)
        disp("Expmt: "+expmt_mode+", Test idx: "+num2str(expmt_id)+", Users: "+num2str(num_users));
        if(expmt_mode=="chan_est")
            %% Collect channels and do consistency checks
            conj_idx = 2;
            zero_chan_idx = mod((conj_idx+1)-1,4)+1;
            chan_ant2_idx_1 = mod((conj_idx+2)-1,4)+1;
            chan_ant2_idx_2 = mod((conj_idx+1)-1,4)+1;
            chan_ant3 = mod((conj_idx+3)-1,4)+1;
            chan_ant5 = mod((conj_idx+2)-1,4)+1;
            chan_ant7 = mod((conj_idx+3)-1,4)+1;
            sup_idx = mod((conj_idx+3)-1,4)+1;
            consistency_indexes = [zero_chan_idx,chan_ant2_idx_1,chan_ant2_idx_2];
            sup_indexes = [chan_ant3,chan_ant5,chan_ant7,sup_idx];
            
            use_pre = false;
            
            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [1,1,1,1];
            if(use_pre)
                [raw_rx_samps_cpo_est, chan_est_cpo_est] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,raw_rx_samps_cpo_est);
            else
                [raw_rx_samps_cpo_est, chan_est_cpo_est] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false);
            end
            
            if(chan_est_cpo_est==-1)
                return
            end
            
            cpo_vec = cpo_est(chan_est_cpo_est,num_users,conj_idx,samp_fac, ofdm_params.SC_IND_DATA);
            
            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [0,0,0,1];
            antenna_state(2,1:4) = [0,1,0,0];
            antenna_state(3,1:4) = [1,0,0,0];
            if(use_pre)
                [raw_rx_samps_A, chan_est_A] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,raw_rx_samps_A);
            else
                [raw_rx_samps_A, chan_est_A] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false);
            end
            
            if(chan_est_A==-1)
                return
            end
            
            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [0,0,0,1];
            antenna_state(4,1:4) = [0,0,1,0];
            antenna_state(5,1:4) = [0,1,0,0];
            antenna_state(6,1:4) = [1,0,0,0];
%             antenna_state(1,1:4) = [0,0,0,1];
%             antenna_state(2,1:4) = [0,0,1,0];
%             antenna_state(3,1:4) = [1,0,0,0];
            
            if(use_pre)
                [raw_rx_samps_B, chan_est_B] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,raw_rx_samps_B);
            else
                [raw_rx_samps_B, chan_est_B] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false);
            end
            
            if(chan_est_B==-1)
                return
            end
            
            
%             antenna_state = zeros([8,4]);
%             antenna_state(1,1:4) = [0,0,0,1];
%             antenna_state(7,1:4) = [0,1,0,0];
%             antenna_state(8,1:4) = [1,0,0,0];
            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [0,0,0,1];
            antenna_state(7,1:4) = [0,0,1,0];
            antenna_state(2,1:4) = [1,0,0,0];
            if(use_pre)
                [raw_rx_samps_C, chan_est_C] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,raw_rx_samps_C);
            else
                [raw_rx_samps_C, chan_est_C] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false);
            end
            
            if(chan_est_C==-1)
                return
            end
            
            
            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [0,0,0,1];
            antenna_state(2,1:4) = [0,0,1,1];
            antenna_state(3,1:4) = [1,1,1,0];
            antenna_state(4,1:4) = [0,1,0,1];
            antenna_state(5,1:4) = [1,0,1,1];
            antenna_state(6,1:4) = [0,1,0,0];
            antenna_state(7,1:4) = [1,0,1,1];
            antenna_state(8,1:4) = [0,1,0,0];
            
            if(use_pre)
                [raw_rx_samps_sup, chan_est_sup] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,raw_rx_samps_sup);
            else
                [raw_rx_samps_sup, chan_est_sup] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false);
            end
            
            if(chan_est_sup==-1)
                return
            end
            
            figure(3)
            plot_chans(chan_est_sup,user_idx, samp_fac, ofdm_params.SC_IND_DATA, conj_idx);
            
            if(disp_plots)
                user_idx=1;
                figure(user_idx)
                consistency_plots(chan_est_A,chan_est_B,chan_est_C,chan_est_sup,user_idx,ofdm_params.SC_IND_DATA,consistency_indexes,cpo_vec,sup_indexes);
            end
            
            rx_samps_all = [raw_rx_samps_A;raw_rx_samps_B;raw_rx_samps_C;raw_rx_samps_sup]; 
%             save(path_dir+"rx_samps_"+rx_samp_file_suffix+".mat",'rx_samps_all');
            
        elseif(expmt_mode=="analog_choice")
        elseif(expmt_mode=="random")
            
        end
    end
end

function cpo_vec = cpo_est(chan_est_cpo_est,num_users,conj_ind,samp_fac,non_nulls)
    cpo_vec = zeros([num_users,samp_fac]);
    cpo_vec_deg = zeros([num_users,samp_fac]);
    disp_plot=true;
    for user_idx=1:1:num_users
        chan_mat_user = squeeze(chan_est_cpo_est(:,user_idx,:,non_nulls));
        ref_phase_mult = exp(-1j*angle(chan_mat_user(:,conj_ind,:)));
        ref_phase_mult_mat = repmat(ref_phase_mult,[1,samp_fac,1]); 
        chan_phasors = exp(1j*angle(ref_phase_mult_mat.*chan_mat_user));
        
        cpo_vec(user_idx,:) = exp(-1j*angle(mean(chan_phasors,[1,3])));
        cpo_vec_deg(user_idx,:) = angle(cpo_vec(user_idx,:))*180/pi;
        
        if(disp_plot)
            for i=1:1:samp_fac
                subplot(2,2,i)
                angles_deg = squeeze(unwrap(angle(chan_phasors(:,i,:)))*180/pi).';
                plot(angles_deg)
                ylim([min(angles_deg,[],'all')-20,max(angles_deg,[],'all')+20])
                hold on
                title("RF Idx: "+num2str(i))
            end
        end
%         cpo_conf(user_idx,:) = std(chan_angles*180/pi,[1,3]);
    end
    
    
end

function plot_chans(chan_est, user_idx, samp_fac, non_nulls, conj_idx)
    chan_est_usr = squeeze(chan_est(:,user_idx,:,:));
    conj_mult = exp(-1j*angle(squeeze(chan_est_usr(:,conj_idx,:)).'));
    for i=1:1:samp_fac
        subplot(4,2,(i-1)*2+1)
        chan_est_rfi = squeeze(chan_est_usr(:,i,:)).';
        conj_chan_est_rfi = chan_est_rfi.*conj_mult;
        chan_db = 20*log10(abs(conj_chan_est_rfi));
        plot(chan_db)
        ylim([min(chan_db(non_nulls,:),[],'all')-10,max(chan_db,[],'all')+10])
        subplot(4,2,(i-1)*2+2)
        plot(unwrap(fftshift(angle(conj_chan_est_rfi(non_nulls,:))))*180/pi)
        ylim([-200,200])
    end
end

function consistency_plots(chan_est_A,chan_est_B,chan_est_C,chan_est_sup,user_idx,non_nulls,consistency_indexes,cpo_vec,sup_indexes)
    % top row: zero consistency
    % middle row: ant1-2 consistency
    % bottom row: sup consistency
%     offset = 0;
%     zero_chan_idx = mod(4-1+offset,4)+1;
    cpo_vec_user = cpo_vec(user_idx,:);
    zero_chan_idx = consistency_indexes(1);
    
    zero_chan_1 = squeeze(chan_est_A(:,user_idx,zero_chan_idx,:)*(cpo_vec_user(consistency_indexes(1)))).';
    zero_chan_2 = squeeze(chan_est_C(:,user_idx,zero_chan_idx+1,:)*(cpo_vec_user(consistency_indexes(1)+1))).';
    subplot(1,2,1)
    plot(20*log10(abs(zero_chan_1)),'Color','b')
    title("Zero chan consistency, mag")
    hold on
    plot(20*log10(abs(zero_chan_2)),'Color','r')
    ylim([-40,20])
    subplot(1,2,2)
    plot(unwrap(fftshift(angle(zero_chan_1(non_nulls,:))))*180/pi,'Color','b')
    title("Zero chan consistency, phase")
    ylim([-200,200])
    hold on
    plot(unwrap(fftshift(angle(zero_chan_2(non_nulls,:))))*180/pi,'Color','r')
    
%     ant2_chan_1 = squeeze(chan_est_A(:,user_idx,consistency_indexes(2),:)*cpo_vec_user(consistency_indexes(2))).';
%     ant2_chan_2 = squeeze(chan_est_sup(:,user_idx,consistency_indexes(3),:)*cpo_vec_user(consistency_indexes(3))).';
%     subplot(3,2,3)
%     plot(20*log10(abs(ant2_chan_1)),'Color','b')
%     title("Ant2 chan consistency, mag")
%     hold on
%     plot(20*log10(abs(ant2_chan_2)),'Color','r')
%     ylim([-40,20])
%     subplot(3,2,4)
%     plot(unwrap(fftshift(angle(ant2_chan_1(non_nulls,:))))*180/pi,'Color','b')
%     title("Ant2 chan consistency, phase")
%     ylim([-200,200])
%     hold on
%     plot(unwrap(fftshift(angle(ant2_chan_2(non_nulls,:))))*180/pi,'Color','r')
%     
%     
%     sup_chan = squeeze(chan_est_sup(1,user_idx,sup_indexes(4),:)*cpo_vec_user(sup_indexes(4))).';
%     sup_chan_A = squeeze(chan_est_A(1,user_idx,sup_indexes(1),:)*cpo_vec_user(sup_indexes(1))).';
%     sup_chan_B = squeeze(chan_est_B(1,user_idx,sup_indexes(2),:)*cpo_vec_user(sup_indexes(2))).';
%     sup_chan_C = squeeze(chan_est_C(1,user_idx,sup_indexes(3),:)*cpo_vec_user(sup_indexes(3))).';
%     zero_chan_sq = squeeze(zero_chan_1(:,1)).';
%     sup_chan_recon = sup_chan_A+sup_chan_B+sup_chan_C-zero_chan_sq;
%     subplot(3,2,5)
%     plot(20*log10(abs(sup_chan)),'Color','b')
%     title("Ant2 chan consistency, mag")
%     hold on
%     plot(20*log10(abs(sup_chan_recon)),'Color','r')
%     ylim([-40,20])
%     subplot(3,2,6)
%     plot(unwrap(fftshift(angle(sup_chan(:,non_nulls))))*180/pi,'Color','b')
%     title("Ant2 chan consistency, phase")
%     ylim([-200,200])
%     hold on
%     plot(unwrap(fftshift(angle(sup_chan_recon(:,non_nulls))))*180/pi,'Color','r')
%     
    
end

