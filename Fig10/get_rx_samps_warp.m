function [rx_samps] = get_rx_samps_warp(antenna_state, bb_gain,rf_gain, chan_tune)
    if(nargin==1)
        bb_gain = 11;
        rf_gain = 3;
        chan_tune = 14;
    elseif(nargin==2)
        rf_gain = 3;
        chan_tune = 14;
        disp("Tuned to channel: "+num2str(chan_tune))
    elseif(nargin==3)
        chan_tune = 14;
        disp("Tuned to channel: "+num2str(chan_tune))
    elseif(nargin==4)
        disp("No default arguments being used, tuning to channel: "+num2str(chan_tune))
    elseif(nargin>4)
        disp("Wrong number of arguments for get_rx_samps_warp function, please check and retry")
        rx_samps=-1;
        return
    end
    %% Configure the PCB switching correctly
    antenna_decimal_states = antenna_state(:,1)*2^3+antenna_state(:,2)*2^2+antenna_state(:,3)*2+antenna_state(:,4);
    antenna_hex_states = fliplr(dec2hex(antenna_decimal_states).');
    all_antenna_dec =  uint32(hex2dec(antenna_hex_states));
    if(all_antenna_dec>0)
        serial_obj = serial('/dev/ttyUSB12','BaudRate',9600);
        fopen(serial_obj);
        pause(0.5);
        fwrite(serial_obj,all_antenna_dec,'uint32');
        pause(0.1)
        rx_str_dec = fread(serial_obj, 1, 'uint32');
        if isempty(rx_str_dec)
            disp("Serial timeout bug :(")
            fclose(serial_obj);
            pause(0.1)
            clear serial_obj;
            return
        end
        rx_str_hex = dec2hex(rx_str_dec)
        if(all_antenna_dec~=rx_str_dec)
            antenna_hex_states
            rx_str_dec
            disp("Antenna state command failed")
            fclose(serial_obj);
            pause(0.1)
            clear serial_obj;
            return
        end
        %% Close serial port
        fclose(serial_obj);
        pause(0.1)
        clear serial_obj;
    end
    %% Now configure the WARP to collect samples
    read_buffs=true;
    save_file=true;
    delayns=6.25*2;

    CHANNEL                 = chan_tune;          % Channel to tune Tx and Rx radios https://warpproject.org/trac/wiki/WARPLab/Reference/Interface/X245
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

        RxGainRF_low = 2;                  % Rx RF Gain in [1:3]
        if(nargin==3)
            disp("Setting gain to: "+num2str(rf_gain)+","+num2str(bb_gain))
        end
        RxGainRF_high = rf_gain;                  % Rx RF Gain in [1:3]
        RxGainBB_low = bb_gain;                 % Rx Baseband Gain in [0:31]
        RxGainBB_high = 22;                 % Rx Baseband Gain in [0:31]
        wl_interfaceCmd(node_rx, ifc_ids_RX.RF_A, 'rx_gains', RxGainRF_low, RxGainBB_high);
        wl_interfaceCmd(node_rx, ifc_ids_RX.RF_B, 'rx_gains', RxGainRF_high, RxGainBB_low);
        wl_interfaceCmd(node_rx, ifc_ids_RX.RF_C, 'rx_gains', RxGainRF_low, RxGainBB_high);
        wl_interfaceCmd(node_rx, ifc_ids_RX.RF_D, 'rx_gains', RxGainRF_high, RxGainBB_low);

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
    %% Collect samples now
    eth_trig.send();
    rx_vec_air_A = wl_basebandCmd(node_rx, RX_RF_VEC_A, 'read_IQ', 0, TX_NUM_SAMPS);
    rx_vec_air_B = wl_basebandCmd(node_rx, RX_RF_VEC_B, 'read_IQ', 0, TX_NUM_SAMPS);
    rx_vec_air_C = wl_basebandCmd(node_rx, RX_RF_VEC_C, 'read_IQ', 0, TX_NUM_SAMPS);
    rx_vec_air_D = wl_basebandCmd(node_rx, RX_RF_VEC_D, 'read_IQ', 0, TX_NUM_SAMPS);

    rx_vec_air_A = rx_vec_air_A(:).';
    rx_vec_air_B = rx_vec_air_B(:).';
    rx_vec_air_C = rx_vec_air_C(:).';
    rx_vec_air_D = rx_vec_air_D(:).';



    rx_samps = [rx_vec_air_A;rx_vec_air_B;rx_vec_air_C;rx_vec_air_D;];
    wl_triggerManagerCmd(nodes, 'output_state_clear', [trig_out_ids.BASEBAND, trig_out_ids.AGC, trig_out_ids.EXT_OUT_P0]);
    wl_basebandCmd(node_rx, RX_RF_ALL, 'tx_rx_buff_dis');
    wl_interfaceCmd(node_rx, RX_RF_ALL, 'tx_rx_dis');
end

