function [lts_ind_first_user] = pack_detect(raw_rx_dec,params,lts_thresh)

unsynch = true;

if(nargin==2)
    lts_thresh=0.5;
end

if(~isfield(params,'max_len_to_check'))
    params.max_len_to_check = min(40000,length(raw_rx_dec));
else
    params.max_len_to_check = min(params.max_len_to_check, length(raw_rx_dec));
end

% Complex cross correlation of Rx waveform with time-domain LTS
lts_corr = abs(conv(conj(fliplr(params.lts_t)), sign(raw_rx_dec(1:params.max_len_to_check))));

% Skip early and late samples - avoids occasional false positives from pre-AGC samples
lts_corr = lts_corr(params.N_SC/2:end-params.N_SC/2);
lts_corr=lts_corr(1:min(params.max_len_to_check,length(lts_corr)));
% Find all correlation peaks
params.LTS_CORR_THRESH=lts_thresh;
% disp("Using lts_thresh: "+num2str(params.LTS_CORR_THRESH))
lts_peaks = find(lts_corr > params.LTS_CORR_THRESH*max(lts_corr));

% Select best candidate correlation peak as LTS-payload boundary
[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);

[lts_last_peak_index,~] = find(((LTS2-LTS1 == length(params.lts_t)*(params.NUM_LTS-1))));

% Set a random value if no valid correlation peak was found
if(isempty(lts_last_peak_index))
%     fprintf('No LTS Correlation Peaks Found!\n');
    lts_ind_first_user = -1;
    return;
end

if(~unsynch)

    sensitivity = 0;
    [LTS1, LTS2] = meshgrid(lts_peaks(lts_last_peak_index),lts_peaks(lts_last_peak_index));
    perfect_condition = (LTS2-LTS1 == (length(params.lts_t)*(params.NUM_LTS+0.5)+params.inter_user_zeros)*(params.num_users-1));
    sensitivity_condition1 = ((LTS2-LTS1) < ((length(params.lts_t)*(params.NUM_LTS+0.5)+params.inter_user_zeros)*(params.num_users-1)+sensitivity));
    sensitivity_condition2 = ((LTS2-LTS1) > ((length(params.lts_t)*(params.NUM_LTS+0.5)+params.inter_user_zeros)*(params.num_users-1)-sensitivity));
    if(sensitivity>0)
        [lts_last_user_index,~] = find(perfect_condition | (sensitivity_condition1 & sensitivity_condition2));
    else
        [lts_last_user_index,~] = find(perfect_condition);
    end



    % Set a random value if no valid correlation peak was found
    if(isempty(lts_last_user_index))
%         fprintf('Multi user LTS detection failed\n');
        lts_ind_first_user = -1;
        return;
    end



    % Set the sample indices of the payload symbols and preamble
    payload_ind = lts_peaks(min(lts_last_peak_index(lts_last_user_index))) + params.N_SC/2+params.inter_user_zeros;
    lts_ind_last_user = payload_ind-(params.NUM_LTS+0.5)*params.N_SC-params.inter_user_zeros;
    lts_ind_first_user = lts_ind_last_user - ((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)*(params.num_users-1);
else
%     lts_range_possibilities = ((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)*(params.num_users)+params.inter_user_zeros;
    lts_range_possibilities_oneuser_min = ((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)-params.inter_user_zeros;
    lts_range_possibilities_oneuser_max = ((params.NUM_LTS+0.5)*params.N_SC+params.inter_user_zeros)+params.inter_user_zeros;
    
    
    sing_user_lts_peaks_abs_idx = lts_peaks(lts_last_peak_index);
    cleaned_peaks =[];
    peak_idx=1;
    while(peak_idx<=numel(sing_user_lts_peaks_abs_idx))
        curr_peak = sing_user_lts_peaks_abs_idx(peak_idx);
        possiblemultipath_peaks = ((sing_user_lts_peaks_abs_idx-curr_peak)>0)...
                & ((sing_user_lts_peaks_abs_idx-curr_peak)<params.N_SC); %Remove peaks within one frame
        if(sum(possiblemultipath_peaks)>0)
%             disp("Multipath peaks found")
        end
        peak_idx=peak_idx+1+sum(possiblemultipath_peaks);
        cleaned_peaks=[cleaned_peaks curr_peak];
    end
        
   sing_user_lts_peaks_abs_idx = cleaned_peaks;
    
    lts_peaks_all_users=[];
    for sing_user_lts_peak_abs_idx=sing_user_lts_peaks_abs_idx
        candidate_peaks = [sing_user_lts_peak_abs_idx];
        peak_find_success=true;
        for other_user_candidates = 1:1:params.num_users-1
            curr_peak = candidate_peaks(end);
            possible_next_indice = ((sing_user_lts_peaks_abs_idx-curr_peak)>=lts_range_possibilities_oneuser_min)...
                & ((sing_user_lts_peaks_abs_idx-curr_peak)<lts_range_possibilities_oneuser_max);
            num_other_possible_peaks = sum(possible_next_indice);
            if(num_other_possible_peaks==1)
                candidate_peaks=[candidate_peaks sing_user_lts_peaks_abs_idx(find(possible_next_indice))];
            elseif(num_other_possible_peaks>1)
%                 disp("Multipeaks found")
                peak_find_success=false;
                break
            elseif(num_other_possible_peaks==0)
                peak_find_success=false;
                break
            end
            
        end
        if(peak_find_success)
            lts_peaks_all_users=[lts_peaks_all_users; candidate_peaks];
        end
    end
    if(isempty(lts_peaks_all_users))
%         disp("Multi user unsynch LTS not found")
        lts_ind_first_user =-1;
        
%         plot(abs(raw_rx_dec))
%         hold on
%         for i=1:1:numel(cleaned_peaks)
%             xline(cleaned_peaks(i))
%         end
        return
    end
%      plot(abs(raw_rx_dec))
%     hold on
%     for i=1:1:4
%         xline(lts_peaks_all_users(1,i))
%     end
%     lts_peaks_all_users-lts_peaks_all_users(:,1)
    lts_ind_first_user = lts_peaks_all_users;
end
% 
% plot_samp_size_r=10000;
% plot_samp_size_l=500;
% plot(abs(raw_rx_dec(lts_ind_first_user-plot_samp_size_l:lts_ind_first_user+plot_samp_size_r)))
% hold on
% xline(plot_samp_size_l)
% xline(payload_ind-(lts_ind_first_user - plot_samp_size_l))

% disabling for now
% cur_last_peak_ind = 2;
% while(lts_ind_last_user<=0)
%     payload_ind = lts_peaks(min(lts_last_peak_index(lts_last_user_index))) + params.N_SC/2+params.inter_user_zeros;
%     lts_ind_last_user = payload_ind-(params.NUM_LTS+0.5)*params.N_SC-params.inter_user_zeros;
%     cur_last_peak_ind = cur_last_peak_ind+1;
% end
% Rx LTS correlation

if(params.plot_flag)
    figure(params.cf);clf;
    lts_to_plot = lts_corr(1:params.max_len_to_check);
    plot(lts_to_plot, '.-b', 'LineWidth', 1);
    hold on;
    grid on;
    line([1 length(lts_to_plot)], params.LTS_CORR_THRESH*max(lts_to_plot)*[1 1], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2);
    title('LTS Correlation and Threshold')
    xlabel('Sample Index')
    myAxis = axis();
    axis([1, 12000, myAxis(3), myAxis(4)]);
end
end
