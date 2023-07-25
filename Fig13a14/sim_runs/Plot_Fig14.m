% pre-computed data
load('./Pre_run_data/nautilus_SimDataFig14.mat')

ant_sweeps = [8,16,32,48,64,128,256];
lwidth=3;
plot_mask = [...% Algo mask description
              1,...% 1: Full Digital
              0,...% 2: FC HBF
              0,...% 3: PC HBF
              0,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              1,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x


% 45,30,20,5
% rfc_pow = 1.2*1000;
bb_pow_const = (((1260)/64)); %4200*0.45
rf_pow_const = (((840)/64)); %4200*0.45
% aa_pow_const_pc = 
cool_pow_tx = 0.5*4200;
fixed_pow_tx = [ ...% Algo mask description
              cool_pow_tx,...% 1: Full Digital
              0,...% 2: FC HBF
              cool_pow_tx+(bb_pow_const+rf_pow_const)*8,...% 3: PC HBF
              0,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              cool_pow_tx+rf_pow_const*8+bb_pow_const*8,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x

cool_pow_rx = 0.05*4200;
fixed_pow_rx = [ ...% Algo mask description
              cool_pow_rx,...% 1: Full Digital
              0,...% 2: FC HBF
              cool_pow_rx+(bb_pow_const+rf_pow_const)*8,...% 3: PC HBF
              0,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              cool_pow_rx+rf_pow_const*8+bb_pow_const*8,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x

var_pow= [ ...% Algo mask description
              bb_pow_const+rf_pow_const,...% 1: Full Digital
              0,...% 2: FC HBF
              0,...% 3: PC HBF
              0,...% 4: UP Dig (DBF 1x)
              0,...% 5: id DBF (Calib)
              0,...% 6: BABF
              0,...% 7: Filed LN
              0,...% 8: DBF 2x
              0,...% 9: DBF 3x
              0];  % 10: DBF 4x



for i=1:1:numel(plot_mask)
    if(plot_mask(i)==1)
        % mean sinrs: evals*sweeps
        
        curr_mean_sinrs_across_users = sinr_mat_evals(:,:,i,1);
        curr_mean_std_across_users = sinr_mat_evals(:,:,i,4);
%         curr_spread_sinrs_across_users = sinr_mat_evals(:,:,i,3)-sinr_mat_evals(:,:,i,2);
        mean_sinr_across_evals = squeeze(mean(curr_mean_sinrs_across_users,1));
%         mean_spread = squeeze(mean(curr_spread_sinrs_across_users,1));
        mean_std_across_evals = squeeze(mean(curr_mean_std_across_users,1));


        figure(3)
        curr_capacity = 100*log2(1+10.^(mean_sinr_across_evals/10))*8;
        net_pow_tx = fixed_pow_tx(i)+var_pow(i).*ant_sweeps;
        net_pow_rx = fixed_pow_rx(i)+var_pow(i).*ant_sweeps;
        EE = curr_capacity./net_pow_tx;
        h = plot(ant_sweeps,EE);
        h.LineStyle="--";
        h.Marker="*";
        h.MarkerSize=10;
        h.LineWidth=lwidth;
        hold on
        
        figure(4)
        SE = curr_capacity/100;
        h = plot(SE,EE);
        h.LineStyle="--";
        h.Marker="*";
        h.MarkerSize=10;
        h.LineWidth=lwidth;
        hold on

        figure(5)
        h = plot(ant_sweeps,SE);
        h.LineStyle="--";
        h.Marker="*";
        h.MarkerSize=10;
        h.LineWidth=lwidth;
        hold on

    end
end


net_str = ["DBF N Ant: N RFC","Fully connected","Partially connected",...
    "DBF 8 RF","N VRF N Ant (Calib, should be same as DBF 1x)","GreenMO N Ant: 8 VRFCs","GreenMO: FilledLN"...
    ,"DBF: 2x","DBF: 3x","DBF: 4x"];
% selected_algos = fliplr(net_str);
selected_algos = net_str(plot_mask==1);

figure(3)
legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("Number of antennas")
ylabel("EE (Mb per joule)")
title("")
grid on
h=xline(64,'-','64 Antenna AAU');
h.LineWidth = 4;
h.HandleVisibility='off';
h.LabelVerticalAlignment='bottom';
h.LabelHorizontalAlignment='center';


figure(4)
legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("SE (bits per hz)")
ylabel("EE (Mb per joule)")
h=xline(50,'-','EE at existing SE level');
h.LineWidth = 4;
h.HandleVisibility='off';
h.LabelVerticalAlignment='middle';
h.LabelHorizontalAlignment='center';
title("")
grid on


figure(5)
legend([selected_algos])
% legend([selected_algos,"SINR requirements"])
xlabel("Number of antennas")
ylabel("SE (Bits per Hz)")
h=yline(50,'-','Current AAU SE level');
h.LineWidth = 4;
h.HandleVisibility='off';
h.LabelHorizontalAlignment='center';
h=xline(64,'-','64 Antenna AAU');
h.LineWidth = 4;
h.HandleVisibility='off';
h.LabelVerticalAlignment='bottom';
h.LabelHorizontalAlignment='center';

title("")
grid on
% figure(2)
% legend([selected_algos])
% % legend([selected_algos,"SINR requirements"])
% xlabel("Number of users")
% ylabel("Energy Efficiency (Mb per Joule)")
% title("")
