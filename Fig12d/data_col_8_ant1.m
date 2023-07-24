root_dir = "/media/wcsng-20/57810C9C5890F0C7/";
savpath_dir = root_dir+"ag_data2/8ant_sep_clk";
dir_prefix = "rx_data";

posn_index = 1;
num_ant = 8;
if(num_ant==8)
    ants_to_use = [1,2,3,4,5,6,7,8];
elseif(num_ant==6)
    ants_to_use = [1,2,3,6,7,8];
elseif(num_ant==4)
    ants_to_use = [1,2,3,4];
end    
max_ants = 8;
expmt_str="rx_data"+num2str(num_ant)+"_"+num2str(posn_index);
% ret_val = system("mkdir "+savpath_dir+"/"+expmt_str);
path_dir = savpath_dir+"/"+expmt_str+"/";
if ~exist(path_dir, 'dir')
   mkdir(path_dir)
end

% expmt_modes = ["analog_choice","random"];
% num_expmts = [1,10];
expmt_modes = ["bit_phase"];
num_expmts = [2,3];
expmt_offset = [0,0];
ch_est_gain = 10;

%% set up ofdm params
ofdm_params = wcsng_ofdm_param_gen(64);
ofdm_params.num_users = 4;
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
disp_plots=false;
samp_fac = 4;


%% final experiment loops
clear ch_est_rx_samps 
for expmt_mode_idx=1:1:numel(expmt_modes)
    expmt_mode_str = expmt_modes(expmt_mode_idx);
    num_trials = num_expmts(expmt_mode_idx);
    curr_expmt_offset = expmt_offset(expmt_mode_idx);
    for expmt_id=(1+curr_expmt_offset):1:num_trials
        disp("Expmt: "+expmt_mode_str+", Test idx: "+num2str(expmt_id)+", Users: "+num2str(num_users));
        if(expmt_mode_str=="bit_phase")
            %% Collect channels and do consistency checks
            conj_idx = 3;

            ant_chan_idx = zeros(1,max_ants);
            ant_chan_num_idx = zeros(1,max_ants);
            
            
            use_pre = false;
            clear rx_samps_ch_est;
            if(use_pre)
                load(path_dir+"rx_samps_ch_est"+num2str(expmt_id)+".mat",'rx_samps_ch_est');
            end
            
            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [1,1,1,1];
            if(use_pre)
                [~, chan_est_cpo_est] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,rx_samps_ch_est{1},ch_est_gain);
            else
                [rx_samps_ch_est{1}, chan_est_cpo_est] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,ch_est_gain);
            end

            if(chan_est_cpo_est==-1)
                disp("Chan est CPO failed, Expmt id: "+num2str(expmt_id))
                return
            end

            cpo_vec = cpo_est(chan_est_cpo_est,num_users,conj_idx,samp_fac, ofdm_params.SC_IND_DATA);

            antenna_state = zeros([8,4]);
            
            antenna_state(1,1:4) = [1,0,0,0];
            antenna_state(2,1:4) = [0,0,1,0];
            antenna_state(3,1:4) = [0,1,0,0];
            ant_chan_idx(1) = conj_idx;
            zero_chan_idx_1 = mod((conj_idx+1)-1,4)+1;
            ant_chan_idx(2) = mod((conj_idx+2)-1,4)+1;
            ant_chan_idx(3) = mod((conj_idx+3)-1,4)+1;
            ant_chan_num_idx(1:3) = ones(1,3);
            if(use_pre)
                [~, chan_est_A] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,rx_samps_ch_est{2},ch_est_gain);
            else
                [rx_samps_ch_est{2}, chan_est_A] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,ch_est_gain);
            end

            if(chan_est_A==-1)
                disp("Chan est A failed, Expmt id: "+num2str(expmt_id))
                return
            end

            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [1,0,0,0];
            antenna_state(4,1:4) = [0,0,0,1];
            antenna_state(5,1:4) = [0,0,1,0];
            antenna_state(6,1:4) = [0,1,0,0];
            ant_chan_idx(4) = mod((conj_idx+1)-1,4)+1;
            ant_chan_idx(5) = mod((conj_idx+2)-1,4)+1;
            ant_chan_idx(6) = mod((conj_idx+3)-1,4)+1;
            ant_chan_num_idx(4:6) = 2*ones(1,3);

            if(use_pre)
                [~, chan_est_B] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,rx_samps_ch_est{3},ch_est_gain);
            else
                [rx_samps_ch_est{3}, chan_est_B] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,ch_est_gain);
            end

            if(chan_est_B==-1)
                disp("Chan est B failed, Expmt id: "+num2str(expmt_id))
                return
            end

            antenna_state = zeros([8,4]);
            antenna_state(1,1:4) = [1,0,0,0];
            antenna_state(7,1:4) = [0,0,0,1];
            antenna_state(8,1:4) = [0,0,1,0];
            ant_chan_idx(7) = mod((conj_idx+1)-1,4)+1;
            ant_chan_idx(8) = mod((conj_idx+2)-1,4)+1;
            zero_chan_idx_2 = mod((conj_idx+3)-1,4)+1;
            ant_chan_num_idx(7:8) = 3*ones(1,2);
            
            if(use_pre)
                [~, chan_est_C] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,rx_samps_ch_est{4},ch_est_gain);
            else
                [rx_samps_ch_est{4}, chan_est_C] = get_chan_est(antenna_state, ofdm_params, samp_fac, conj_idx, num_lts_to_use,false,ch_est_gain);
            end

            if(chan_est_C==-1)
                disp("Chan est C failed, Expmt id: "+num2str(expmt_id))
                return
            end

            if(disp_plots)
                consistency_indexes=[zero_chan_idx_1,zero_chan_idx_2];
                for user_idx=1:1:num_users
                    figure(user_idx)
                    consistency_plot(chan_est_A,chan_est_C,user_idx,ofdm_params.SC_IND_DATA,consistency_indexes,cpo_vec);
                end
            end
            
            num_subc =ofdm_params.N_SC;
            cpo_vec_mat = repmat(cpo_vec,[1,1,num_subc]);
            clear chan_est_dict
            chan_est_dict{1} = squeeze(chan_est_A(end,:,:,:)).*cpo_vec_mat;
            chan_est_dict{2} = squeeze(chan_est_B(end,:,:,:)).*cpo_vec_mat;
            chan_est_dict{3} = squeeze(chan_est_C(end,:,:,:)).*cpo_vec_mat;
            
            channel_mat = zeros(num_users,max_ants,num_subc);
            zero_chan = squeeze(chan_est_A(end,:,zero_chan_idx_1,:));
            for i=1:1:max_ants
                curr_chan_idx = ant_chan_num_idx(i);
                curr_pos_idx = ant_chan_idx(i);
                curr_chan_obj = chan_est_dict{curr_chan_idx};
                channel_mat(:,i,:) = squeeze(curr_chan_obj(:,curr_pos_idx,:))-zero_chan;
            end
%             channel_mat
            narrow_beamformer_chan = squeeze(channel_mat(:,:,ofdm_params.SC_IND_DATA(1)));
%             beamnull_mat_all_users = get_beamnull_mat(narrow_beamformer_chan,num_ant);
            beamform_mat_all_users = narrow_beamformer_chan'; % conj+transpose
            
            snr_corr_phase = snr_calc(channel_mat,narrow_beamformer_chan);
%             snr_corr_phase_null = snr_calc(channel_mat,beamnull_mat_all_users);
            
            beamf_antenna_state = analog_beamf_config_choser(beamform_mat_all_users,num_ant,ants_to_use);
            snr_bin_phase = snr_calc(channel_mat,beamf_antenna_state.');
            
            raw_rx_samps = get_rx_samps_warp(beamf_antenna_state,8);
            
%             beamf_antenna_state = analog_beamf_config_choser(beamnull_mat_all_users.',num_ant);
%             snr_bin_phase_null = snr_calc(channel_mat,beamf_antenna_state.');
            
%             raw_rx_samps_beamnull = get_rx_samps_warp(beamf_antenna_state);
            
            save(path_dir+"rx_samps_ch_est_"+num2str(expmt_id)+".mat",'rx_samps_ch_est');
            save(path_dir+"rx_samps_"+expmt_mode_str+"_"+num2str(expmt_id)+".mat",'raw_rx_samps');
%             save(path_dir+"rx_samps_"+expmt_mode_str+"_null"+"_"+num2str(expmt_id)+".mat",'raw_rx_samps_null');
            save(path_dir+"antenna_state_"+expmt_mode_str+"_"+num2str(expmt_id)+".mat",'beamf_antenna_state');
%       
        elseif(expmt_mode_str=="random")
            non_invert_chan_mat = true;
            rand_antenna_state = zeros([8,samp_fac]);
            while(non_invert_chan_mat)
                antenna_state_numant = randi(2,num_ant,samp_fac)-1;
                if(rank(antenna_state_numant)==samp_fac)
                    non_invert_chan_mat=false;
                    break
                end
            end
            % collect samples, process later
            rand_antenna_state(ants_to_use,:) = antenna_state_numant;
            raw_rx_samps = get_rx_samps_warp(rand_antenna_state);
            save(path_dir+"rx_samps_"+expmt_mode_str+"_"+num2str(expmt_id)+".mat",'raw_rx_samps');
            save(path_dir+"antenna_state_"+expmt_mode_str+"_"+num2str(expmt_id)+".mat",'rand_antenna_state');
        end
    end
end

function beamnull_vec = get_beamnull_mat(chan_mat,num_ants_to_use)
    [num_users,max_ants] = size(chan_mat);
    beamnull_vec = zeros(num_users,max_ants);
    chan_mat_to_use = chan_mat(:,1:1:num_ants_to_use);
    for curr_user_idx = 1:1:num_users
        interferers = 1:1:num_users;
        interferers(curr_user_idx)=[];
        chan_to_cancel = chan_mat_to_use(interferers,:);
        null_chan = null(chan_to_cancel);
        sig_chan = chan_mat_to_use(curr_user_idx,:);
        null_comb_vec = sig_chan*null_chan;
        fin_comb_vec = null_chan*null_comb_vec';
        beamnull_vec(curr_user_idx,1:1:num_ants_to_use) = fin_comb_vec;
    end
end

function switched_beamf_mat = analog_beamf_config_choser(beamf_mat_all_users,num_ants_to_use,ants_to_use)
    [max_ants,num_users] = size(beamf_mat_all_users);
    mat_to_use = beamf_mat_all_users(ants_to_use,:);
    non_invert_chan_mat = true;
    ant_idx_perf=1;
    perf_choice = ones(1,num_users);
    while(non_invert_chan_mat)
        switched_beamf_mat= zeros(max_ants,num_users);
        if(sum(perf_choice>2)>0)
            break
        end
        for user_idx=1:1:num_users
            antenna_score_func = zeros(1,num_ants_to_use);
            antenna_configs = zeros(num_ants_to_use,num_ants_to_use);
            curr_vec = mat_to_use(:,user_idx);
            for ant_idx=1:1:num_ants_to_use
                curr_phase = mat_to_use(ant_idx,user_idx);
                conj_vec = curr_vec*conj(curr_phase);
                antenna_configs(ant_idx,:) = abs(angle(conj_vec))<pi/4;
                antenna_score_func(ant_idx) = -sum(antenna_configs(ant_idx,:)); % sort does asc
            end
            [~,antenna_perf_order] = sort(antenna_score_func);
            switched_beamf_mat(ants_to_use,user_idx) = antenna_configs(antenna_perf_order(perf_choice(user_idx)),:);
        end
        
        rank_val = rank(switched_beamf_mat);
        if(rank_val==num_users)
            non_invert_chan_mat=false;
            break
        else
            perf_choice(ant_idx_perf)=perf_choice(ant_idx_perf)+1;
            ant_idx_perf=mod((ant_idx_perf+1)-1,num_users)+1;
        end
            
    end
    
    
    while(non_invert_chan_mat)
        switched_beamf_mat(ants_to_use,:) = randi(2,num_ants_to_use,num_users)-1;
        if(rank(switched_beamf_mat)==num_users)
            non_invert_chan_mat=false;
            disp("Unable to find invertible analog soln, supplying random mat")
            break
        end
    end

    
end

function snr_est = snr_calc(chan_mat, beamf_vec)
    snr_est = zeros(4,4);
    narrow_est = zeros(4,4);
%     figure(12)
    for user_idx=1:1:4
        curr_beamf_vec = conj(beamf_vec(user_idx,:));
        curr_beamf_mat = repmat(curr_beamf_vec,[4,1,64]);
        chan_mult = chan_mat.*curr_beamf_mat;
        snr_est(user_idx,:) = sum(abs(sum(chan_mult,[2])),[3]);
        subc_profile = abs(sum(chan_mult,[2]));
        narrow_est(user_idx,:) = subc_profile(:,2);
%         subplot(2,2,user_idx)
%         plot(
    end
end
function cpo_vec = cpo_est(chan_est_cpo_est,num_users,conj_ind,samp_fac,non_nulls)
    cpo_vec = zeros([num_users,samp_fac]);
    cpo_vec_deg = zeros([num_users,samp_fac]);
    disp_plot=false;
    for user_idx=1:1:num_users
        if(disp_plot)
            figure(user_idx)
        end
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

function consistency_plot(chan_est_A,chan_est_C,user_idx,non_nulls,consistency_indexes,cpo_vec)
    % top row: zero consistency
    % middle row: ant1-2 consistency
    % bottom row: sup consistency
%     offset = 0;
%     zero_chan_idx = mod(4-1+offset,4)+1;
    cpo_vec_user = cpo_vec(user_idx,:);
    zero_chan_idx_1 = consistency_indexes(1);
    zero_chan_idx_2 = consistency_indexes(2);
    
    zero_chan_1 = squeeze(chan_est_A(:,user_idx,zero_chan_idx_1,:)*(cpo_vec_user(zero_chan_idx_1))).';
    zero_chan_2 = squeeze(chan_est_C(:,user_idx,zero_chan_idx_2,:)*(cpo_vec_user(zero_chan_idx_2))).';
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
   
end

