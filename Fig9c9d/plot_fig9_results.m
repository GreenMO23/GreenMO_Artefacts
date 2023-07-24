% The Script takes about 2 minutes to run
clear all
data_dir = "C:\ag_data2\trace_evals";
config = "dbf"; % for Fig 9c
%config = "hbf"; % for Fig 9d
if(config == "dbf")
    expmt_inds = [1,2,3,0];
    legend_str = ["DBF: 4 ant","DBF: 6 ant","DBF: 8 ant","GreenMO"];
else
    expmt_inds = [4,5,6,0];
    legend_str = ["HBF 8-2-4","HBF 8-4-4","HBF 8-8-4","GreenMO"];
end

% Evaluated at 10 user configs
posn_indexes = 1:1:10;



% Loop through experiments: BABF, DBF or HBF depending on
% the config setup above
for expmt_idx = 1:1:numel(expmt_inds)
    % Store SINRs here
    all_sinrs = [];
    % Loop through the 10 positions
    expmt_ind = expmt_inds(expmt_idx);
    disp("Beginning to process Data for: "+expmt_mode)
    for posn_index = posn_indexes
        % Each config is tested 3 times
        test_ids=1:1:3;
        
        if(expmt_ind ==0)
            avg_sinr_across_users = analyse_data_greenmo(posn_index, data_dir, test_ids);
        else
           avg_sinr_across_users = analyse_data_trace_level(posn_index, data_dir, test_ids,expmt_ind); 
        end
        

        all_sinrs = [all_sinrs avg_sinr_across_users];
        expmt_mode = legend_str(expmt_idx);
    
    end
    disp("Finished Processing Data for: "+expmt_mode)
    
    figure(1)
    ph=cdfplot(all_sinrs);
    ph.LineWidth = 4;
    hold on
end
legend(legend_str);
ylabel("CDF");
xlabel("SINRs (dB)");