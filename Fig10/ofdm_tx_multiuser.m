clear ofdm_params_tx
rng(1)
root_dir = "/media/wcsng-20/57810C9C5890F0C7/";
savpath_dir = root_dir+"conf_room_fin/8ant_multipos_QAM16_34r";
num_users = 1;
% expmt_str="twouser_test";
path_dir = savpath_dir+"/tx_data"+num2str(num_users)+"/";
if ~exist(path_dir, 'dir')
       mkdir(path_dir)
end


ofdm_params = wcsng_ofdm_param_gen(64);
ofdm_params.num_packets = 400;
ofdm_params.packet_step_size = 1;
ofdm_params.MOD_ORDER = 16;
ofdm_params.channel_coding_rate = 0.75;
ofdm_params.NUM_LTS = 10;
ofdm_params.num_users = num_users;
ofdm_params.enable_channel_codes=true;

figure(1)

for i=1:1:ofdm_params.num_users
    ofdm_params.user_index = i-1;
    filename_tx = path_dir+"uhd_ofdm_tx_"+num2str(ofdm_params.num_users)+"user_coded_QAM16_34r_"+num2str(i-1)+".dat";
    [tx_samples,ofdm_params_tx(i)] = ofdm_tx(ofdm_params); % Take tx packet from OFDM_tx.m
    repmat_n=1000;
    plot(abs(tx_samples))
%     rms(tx_samples)
    hold on
%     tx_packet_tosave = repmat(tx_samples,repmat_n,1);
%     size(tx_packet_tosave)
%     write_complex_binary(tx_packet_tosave,filename_tx);
%     clear tx_packet_tosave
end
% save(path_dir+"ofdm_params_tx.mat",'ofdm_params_tx');