clear all
% 10 testing configurations
posn_arr = [1,2,3,4,5,6,7,8,9,10];
% Sweep antennas from 4 to 8
num_ants = [8,6,4];
% Base Folder
data_dir = "C:\ag_data2\goodput_evals";
% data_dir = "C:\Users\agrim\Downloads\ag_data2\8ant_multipos";
% 3 tests perconfig 
test_ids=1:1:3;

num_expmts = numel(test_ids);
num_tests = numel(num_ants);
% metrics needed
avg_evm_snr_users_across_tests = zeros(num_tests,num_expmts,4);
avg_sinr_users_across_tests = zeros(num_tests,num_expmts,4);
avg_ber_users_across_tests = zeros(num_tests,num_expmts,4);

% sweep posn
for posn_index = posn_arr
    % sweep test id
    for num_ant = num_ants
        

        % store metrics for later processing (in decode_out)
        % this function stores the data into separate file for later
        % processing
        [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
        avg_sinr_users_across_packets,min_sinr_users_across_packets,...
        avg_ber_users_across_packets,max_ber_users_across_packets]...
        = analyse_data_greenmo(posn_index,num_ant, data_dir, test_ids);
    
    end
    disp("Processed Data for, posn_index "+num2str(posn_index)+", num ants: "+num2str(num_ant))
end

