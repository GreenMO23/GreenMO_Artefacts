posn_index = 9;

folder_name = "8ant_multipos_QAM16_34r";
test_ids=1:1:3;
num_expmts = numel(test_ids);
% num_ants = [1,8,6];
% expmt_modes =["fdma","bit_phase"];
num_ants = [1,8,6,4,4];
expmt_modes =["fdma","bit_phase","bit_phase","bit_phase","identity"];
% num_ants = [1,8,4];
% expmt_modes =["fdma","bit_phase","identity"];
% num_ants = [1,4,4];
% expmt_modes =["fdma","bit_phase","identity"];
num_tests = numel(num_ants);

fdma_ind = 1;

avg_evm_snr_users_across_tests = zeros(num_tests,num_expmts,4);
avg_sinr_users_across_tests = zeros(num_tests,num_expmts,4);
avg_ber_users_across_tests = zeros(num_tests,num_expmts,4);
if(posn_index==9)
    lts_thresh=0.5;
else
    lts_thresh=0.6;
end

for test_ind=1:1:num_tests
    expmt_mode = expmt_modes(test_ind);
    num_ant = num_ants(test_ind);
    if(expmt_mode~="fdma")
        [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
        avg_sinr_users_across_packets,min_sinr_users_across_packets,...
        avg_ber_users_across_packets,max_ber_users_across_packets]...
        = analyse_data_osamp(posn_index,num_ant,expmt_mode, folder_name, test_ids,lts_thresh);
    else
        [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
        avg_sinr_users_across_packets,min_sinr_users_across_packets,...
        avg_ber_users_across_packets,max_ber_users_across_packets]...
        = analyse_data_fdma(posn_index,folder_name, test_ids);
    end

    avg_evm_snr_users_across_tests(test_ind,:,:) = avg_evm_snr_users_across_packets;
    avg_sinr_users_across_tests(test_ind,:,:) = avg_sinr_users_across_packets;
    avg_ber_users_across_tests(test_ind,:,:) = avg_ber_users_across_packets;
    disp("Processed Data for: "+expmt_mode+", num ants: "+num2str(num_ant))
%     avg_evm_snr_users_across_packets
%     min_evm_snr_users_across_packets
    avg_sinr_users_across_packets
%     min_sinr_users_across_packets
%     avg_ber_users_across_packets
%     max_ber_users_across_packets
end
test_inds = 1:1:num_tests;
% test_inds(fdma_ind)=[];
capacity_others  = sum(squeeze(mean(avg_sinr_users_across_tests(test_inds,:,:),2)),2);
capacity_fdma  = sum(mean(avg_sinr_users_across_tests(fdma_ind,:,:),2));
capacity_gain_sinr = (capacity_others/capacity_fdma)

% avg_evm_snrs = squeeze(mean(avg_evm_snr_users_across_tests,2))
avg_ber = squeeze(mean(avg_ber_users_across_tests,2))
