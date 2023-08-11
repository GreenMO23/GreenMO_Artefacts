% Takes about 2 minutes
clear all;
data_dir = "C:\ag_data2\capacity_result";
% 10 positions
posn_arr = 1:1:10;
% 3 tests each
test_ids=1:1:3;
num_expmts = numel(test_ids);
% fdma is single antenna
num_ants = [1,8,6,4];
num_tests = numel(num_ants);
fdma_ind = 1;
all_capacities=zeros(numel(posn_arr),num_tests);

% loop through positions
for posn_index=posn_arr
    
    disp("Processing for position index: "+num2str(posn_index))
    avg_evm_snr_users_across_tests = zeros(num_tests,num_expmts,4);
    avg_sinr_users_across_tests = zeros(num_tests,num_expmts,4);
    avg_ber_users_across_tests = zeros(num_tests,num_expmts,4);

    for test_ind=1:1:num_tests
        num_ant = num_ants(test_ind);
        if(num_ant~=1) % greenmo is multi antenna
            [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
            avg_sinr_users_across_packets,min_sinr_users_across_packets,...
            avg_ber_users_across_packets,max_ber_users_across_packets]...
            = analyse_data_greenmo(posn_index,num_ant,data_dir, test_ids);
        else % fdma is single antenna
            [avg_evm_snr_users_across_packets,min_evm_snr_users_across_packets,...
            avg_sinr_users_across_packets,min_sinr_users_across_packets,...
            avg_ber_users_across_packets,max_ber_users_across_packets]...
            = analyse_data_fdma(posn_index,data_dir, test_ids);
        end


        avg_sinr_users_across_tests(test_ind,:,:) = 10.^(avg_sinr_users_across_packets/10);

    end
    test_inds = 1:1:num_tests;
    
    BW = 10^7;

    all_capacities(posn_index,:)  = BW*sum(squeeze(log2(1+mean(avg_sinr_users_across_tests,2))),2);

end
%% Plot
avg_capacities_not_norm = mean(all_capacities,1);
std_capacities_not_norm = std(all_capacities,0,1);
osamp_caps = all_capacities(:,2:1:4);
fdma_cap = repmat(all_capacities(:,1),[1,3]);
norm_capacities = osamp_caps./fdma_cap;
avg_capacities = squeeze(mean(norm_capacities,1));
std_capacities = squeeze(std(norm_capacities,0,1));

linew=3;
figure(1)
ant_vec = [8,6,4];
ph=errorbar(ant_vec, avg_capacities, std_capacities,'Marker','x','MarkerSize',12);
ph.LineWidth = linew;
ph.Color = "blue";
hold on
ph = yline(1,'-',"FDMA Capacity: 220\pm 33 Mbps",'LabelHorizontalAlignment','left');
ph.LineWidth = linew;
ph.Color = "black";

ylabel("Capacity (Mbps)")
xlabel("Number of Antennas")

ylim([0,1.4])
xlim([3,9])
yticks(0:0.1:1.4)
grid minor
grid on