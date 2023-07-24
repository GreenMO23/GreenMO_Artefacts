%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wl_example_siso_ofdm_txrx.m
% A detailed write-up of this example is available on the wiki:
% http://warpproject.org/trac/wiki/WARPLab/Examples/OFDM
%
% Copyright (c) 2015 Mango Communications - All Rights Reserved
% Distributed under the WARP License (http://warpproject.org/license)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dir = "/media/wcsng-20/57810C9C5890F0C7/";
savpath_dir = root_dir+"ag_data/MMSE_conf_room_QAM16_12r";
num_users = 2;
posn_index = 9;
expmt_str="rx_data"+num2str(num_users)+"_"+num2str(posn_index);
ret_val = system("mkdir "+savpath_dir+"/"+expmt_str);
path_dir = savpath_dir+"/"+expmt_str+"/";
expmt_modes = ["pcb","mmse"];
expmt_ids=[0,1,2];

for expmt_mode=expmt_modes
    serial_obj = serial('/dev/ttyUSB11','BaudRate',9600)
    fopen(serial_obj);
    pause(1);

    if(expmt_mode=="pcb")
        fprintf(serial_obj,'%s\n',"1");
        % fprintf(s,"1\n");
        % Read the applied force
        ser_msg = fgets(serial_obj);
        disp(ser_msg)

        if(ser_msg(1:3)=="PCB")
            disp("PCB config command started successfully, continuing")
            pause(1)
        end

    elseif(expmt_mode=="mmse")
        fprintf(serial_obj,'%s\n',"2");
        % fprintf(s,"2\n");
        % Read the applied force
        ser_msg = fgets(serial_obj);
        disp(ser_msg)

        if(ser_msg(1:4)=="MMSE")
            disp("MMSE command started succesfully, continuing")
            pause(1)
        end
    else
        disp("Invalid mode, exiting code")
        return
    end

    fclose(serial_obj);
    clear serial_obj;

    %%
    % expmt_ids=[0];
    read_buffs=true;
    save_file=true;
    delayns=0;

    CHANNEL                 = 12;          % Channel to tune Tx and Rx radios https://warpproject.org/trac/wiki/WARPLab/Reference/Interface/X245
    % WARPLab experiment params
    USE_AGC                 = false;        % Use the AGC if running on WARP hardware
    MAX_TX_LEN              = 2^20;        % Maximum number of samples to use for this experiment
    TRIGGER_OFFSET_TOL_NS   = 3000;        % Trigger time offset toleration between Tx and Rx that can be accomodated

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the WARPLab experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NUMNODES = 1;

    % Create a vector of node objects
    nodes   = wl_initNodes(NUMNODES);
    node_rx = nodes(1);

    % Create a UDP broadcast trigger and tell each node to be ready for it
    eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(nodes, 'add_ethernet_trigger', [eth_trig]);

    % Read Trigger IDs into workspace
    trig_in_ids  = wl_getTriggerInputIDs(nodes(1));
    trig_out_ids = wl_getTriggerOutputIDs(nodes(1));

    % For both nodes, we will allow Ethernet to trigger the buffer baseband and the AGC
    % wl_triggerManagerCmd(nodes, 'output_config_input_selection', [trig_out_ids.BASEBAND, trig_out_ids.AGC], [trig_in_ids.ETH_A]);
    % wl_triggerManagerCmd(nodes, 'output_config_input_selection', [trig_out_ids.BASEBAND, trig_out_ids.AGC, trig_out_ids.EXT_OUT_P0], [trig_in_ids.ETH_A]);

    % Let ethernet trigger EXT_OUT_P0
    wl_triggerManagerCmd(nodes, 'output_config_input_selection', [trig_out_ids.EXT_OUT_P0], [trig_in_ids.ETH_A]);
    % Wait for trigger back from PCB
    wl_triggerManagerCmd(nodes, 'output_config_input_selection', [trig_out_ids.BASEBAND, trig_out_ids.AGC], [trig_in_ids.EXT_IN_P0]);
    % wl_triggerManagerCmd(nodes, 'input_config_debounce_mode', [trig_in_ids.EXT_IN_P0], 'enable'); 

    wl_triggerManagerCmd(nodes, 'output_config_hold_mode', [trig_out_ids.EXT_OUT_P0], 'enable');

    % Set the trigger output delays.
    nodes.wl_triggerManagerCmd('input_config_delay', [trig_in_ids.EXT_IN_P0], delayns);
    %     nodes.wl_triggerManagerCmd('output_config_delay', [trig_out_ids.AGC], TRIGGER_OFFSET_TOL_NS);

    % Get IDs for the interfaces on the boards. 
    ifc_ids_RX = wl_getInterfaceIDs(node_rx);

    % Set up the TX / RX nodes and RF interfaces
     % receive on RF_A
    RX_RF_A     = ifc_ids_RX.RF_A;
    RX_RF_VEC_A = ifc_ids_RX.RF_A;
    RX_RF_B     = ifc_ids_RX.RF_B;
    RX_RF_VEC_B = ifc_ids_RX.RF_B;
    RX_RF_C     = ifc_ids_RX.RF_C;
    RX_RF_VEC_C = ifc_ids_RX.RF_C;
    RX_RF_D     = ifc_ids_RX.RF_D;
    RX_RF_VEC_D = ifc_ids_RX.RF_D;
    RX_RF_ALL = ifc_ids_RX.RF_ALL;

    % Set up the interface for the experiment
    % 2.4 here means selecting 2.4 GHz
    wl_interfaceCmd(node_rx, RX_RF_ALL, 'channel', 2.4, CHANNEL);

    if(USE_AGC)
        disp("Using AGC")
        wl_interfaceCmd(node_rx, RX_RF_ALL, 'rx_gain_mode', 'automatic');
        wl_basebandCmd(nodes, 'agc_target', -13);
    else
         if(expmt_mode=="pcb")
            RxGainRF_PCB = 3;                  % Rx RF Gain in [1:3]
            RxGainRF = 2;                  % Rx RF Gain in [1:3]
%              RxGainBB = 11;                 % Rx Baseband Gain in [0:31]
            RxGainBB = 10+3;                 % Rx Baseband Gain in [0:31]
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_A, 'rx_gains', RxGainRF, RxGainBB);
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_B, 'rx_gains', RxGainRF_PCB, RxGainBB);
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_C, 'rx_gains', RxGainRF, RxGainBB);
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_D, 'rx_gains', RxGainRF, RxGainBB);
         else
            RxGainRF_1 = 3;                  % Rx RF Gain in [1:3]
            RxGainRF_2 = 3;                  % Rx RF Gain in [1:3]
%             RxGainBB_1 = 26;                 % Rx Baseband Gain in [0:31]
%             RxGainBB_2 = 17;
%             RxGainBB_3 = 15;
            RxGainBB_2 = 10+3;
            RxGainBB_3 = 10+3;
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_A, 'rx_gains', RxGainRF_2, RxGainBB_2);
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_B, 'rx_gains', RxGainRF_1, RxGainBB_3);
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_C, 'rx_gains', RxGainRF_2, RxGainBB_2);
            wl_interfaceCmd(node_rx, ifc_ids_RX.RF_D, 'rx_gains', RxGainRF_1, RxGainBB_3);
         end

    end

    wl_interfaceCmd(node_rx, RX_RF_ALL, 'rx_lpf_corn_freq', 3);
    wl_interfaceCmd(node_rx, RX_RF_ALL, 'rx_lpf_corn_freq_fine', 5);

    % Get parameters from the node
    %     SAMP_FREQ    = wl_basebandCmd(nodes(1), 'tx_buff_clk_freq');
    SAMP_FREQ = 40e6;
    Ts           = 1/SAMP_FREQ;

    example_mode_string = 'hw';



    wl_interfaceCmd(node_rx, RX_RF_A, 'rx_en');
    wl_interfaceCmd(node_rx, RX_RF_B, 'rx_en');
    wl_interfaceCmd(node_rx, RX_RF_C, 'rx_en');
    wl_interfaceCmd(node_rx, RX_RF_D, 'rx_en');
    wl_basebandCmd(node_rx, RX_RF_A, 'rx_buff_en');
    wl_basebandCmd(node_rx, RX_RF_B, 'rx_buff_en');
    wl_basebandCmd(node_rx, RX_RF_C, 'rx_buff_en');
    wl_basebandCmd(node_rx, RX_RF_D, 'rx_buff_en');

    % Retrieve the received waveform from the Rx node
    TX_NUM_SAMPS = 32768*16;
    wl_basebandCmd(node_rx, RX_RF_VEC_A, 'rx_length', TX_NUM_SAMPS);
    wl_basebandCmd(node_rx, RX_RF_VEC_B, 'rx_length', TX_NUM_SAMPS);
    wl_basebandCmd(node_rx, RX_RF_VEC_C, 'rx_length', TX_NUM_SAMPS);
    wl_basebandCmd(node_rx, RX_RF_VEC_D, 'rx_length', TX_NUM_SAMPS);



    for expmt_id=expmt_ids
        disp("Starting Expmt: "+num2str(expmt_id))

        % Trigger the Tx/Rx cycle at both nodes
        eth_trig.send();
        % pause(0.001)
        % wl_triggerManagerCmd(nodes, 'output_state_clear', [trig_out_ids.BASEBAND, trig_out_ids.AGC, trig_out_ids.EXT_OUT_P0]);
        %%
        trig_debug=false;
        if(trig_debug)
            for i=1:1:20
                eth_trig.send();
                pause(0.01)
            end 

        end

        %%
    %     syms_tx = TX_NUM_SAMPS + (ceil((TRIGGER_OFFSET_TOL_NS*1e-9) / (1/SAMP_FREQ)));
        if(read_buffs)
            rx_vec_air_A = wl_basebandCmd(node_rx, RX_RF_VEC_A, 'read_IQ', 0, TX_NUM_SAMPS);
            rx_vec_air_B = wl_basebandCmd(node_rx, RX_RF_VEC_B, 'read_IQ', 0, TX_NUM_SAMPS);
            rx_vec_air_C = wl_basebandCmd(node_rx, RX_RF_VEC_C, 'read_IQ', 0, TX_NUM_SAMPS);
            rx_vec_air_D = wl_basebandCmd(node_rx, RX_RF_VEC_D, 'read_IQ', 0, TX_NUM_SAMPS);

            rx_vec_air_A = rx_vec_air_A(:).';
            rx_vec_air_B = rx_vec_air_B(:).';
            rx_vec_air_C = rx_vec_air_C(:).';
            rx_vec_air_D = rx_vec_air_D(:).';



            rx_samps = [rx_vec_air_A;rx_vec_air_B;rx_vec_air_C;rx_vec_air_D;];
            if(save_file)
                save(path_dir+expmt_mode+"_"+num2str(expmt_id)+".mat",'rx_samps')
            end
        end
        pause(0.5)
        wl_triggerManagerCmd(nodes, 'output_state_clear', [trig_out_ids.BASEBAND, trig_out_ids.AGC, trig_out_ids.EXT_OUT_P0]);
        pause(0.5)
    end
    %%
    wl_triggerManagerCmd(nodes, 'output_state_clear', [trig_out_ids.BASEBAND, trig_out_ids.AGC, trig_out_ids.EXT_OUT_P0]);
    wl_basebandCmd(node_rx, RX_RF_ALL, 'tx_rx_buff_dis');
    wl_interfaceCmd(node_rx, RX_RF_ALL, 'tx_rx_dis');
    pause(0.5)
end

