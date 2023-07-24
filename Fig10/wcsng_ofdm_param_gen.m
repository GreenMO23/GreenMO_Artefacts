function [params] = wcsng_ofdm_param_gen(N_SC)
%OFDM_PARAM_GEN_V2 Generates default parameters for the wcsng ofdm
%library function
%   
if(nargin<1)
    N_SC = 64;
end

    % General Params:
    params.WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
    params.plot_flag               = 0;           % Do you want plots?
    params.CHANNEL                 = 11;          % Channel to tune Tx and Rx radios
    params.calc_stats              = 1;
    params.cf                      = 201;         % Figure Number
    
    % Waveform params
    params.N_OFDM_SYMS             = 50;         % Number of OFDM symbols
    params.MOD_ORDER               = 4;           % Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)
    params.TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1])
    params.N_ZERO_PAD              = 2000;
    params.NUM_LTS                 = 2;
    
    % OFDM params
    [params] = ofdm_param_gen_arb(N_SC,params);
    

    % Rx processing params
    params.FFT_OFFSET                    = 4;           % Number of CP samples to use in FFT (on average)
    params.LTS_CORR_THRESH               = 0.8;         % Normalized threshold for LTS correlation
    params.DO_APPLY_CFO_CORRECTION       = 1;           % Enable CFO estimation/correction
    params.DO_APPLY_PHASE_ERR_CORRECTION = 1;           % Enable Residual CFO estimation/correction
    params.DO_APPLY_SFO_CORRECTION       = 1;           % Enable SFO estimation/correction
    params.DO_DECODE                     = 1;           % Whether to Decode Data or Not

    % Use sane defaults for hardware-dependent params in sim-only version
    params.example_mode_string = 'sim';

    % Channel Coding
    params.enable_channel_codes = false;
    params.channel_coding_rate = .5; % coding rate 
    params.trellis_end_length = 8; % bits for trellis to end 
    
    % MIMO Parameters
    params.mimo_tx = 2; % Max 4 or 5 don't know what SPATIAL STREAM SHIFT DOES
    params.mimo_rx = 4;
    params.num_strm = 1; % Number of TX Streams, strm =1 -> tx diversity
    params.TX_SPATIAL_STREAM_SHIFT = 3;
    
    params.num_strm_rx = 1; % Combining Strategy
    
    % RX only params
    params.num_packets = 1; % Decode only one packet
    params.packet_step_size = 1; % Number of packets to skip between successive decoding
    
    
    params.num_users=1;
    params.user_index=0;
    params.inter_user_zeros=512;
   
    
end

function [params] = ofdm_param_gen_arb(N_SC,params)

    % OFDM params
    params.N_SC                    = N_SC;                                     % Number of subcarriers
    params.SAMP_FREQ               = 20e6;
    params.CP_LEN                  = round(16/20e6*params.SAMP_FREQ);        % Cyclic prefix length
    params.N_STS                   = 10;                                     % Number of STS
    
    subc_spacing = params.SAMP_FREQ/params.N_SC;
    num_filled_sc = round(params.N_SC*52/64);
    num_filled_sc = num_filled_sc + mod(num_filled_sc,2);

    num_pilots_needed = round(num_filled_sc/13);
    num_pilots_needed = num_pilots_needed + mod(num_pilots_needed,2);
    
    FILLED_SC_IND = [2:1:(num_filled_sc/2+1) (params.N_SC-num_filled_sc/2+1):params.N_SC];
    params.FILLED_SC_IND = FILLED_SC_IND;
    
    left_cut = linspace(2,(num_filled_sc/2+1),(num_pilots_needed/2)+1);
    right_cut = linspace((params.N_SC-num_filled_sc/2+1),params.N_SC,(num_pilots_needed/2)+1);
    left_mean = round((left_cut(1:end-1)+left_cut(2:end))/2);
    right_mean = round((right_cut(1:end-1)+right_cut(2:end))/2);
    
    params.SC_IND_PILOTS           = [left_mean right_mean];                           % Pilot subcarrier indices
    pilot_screw = [1,1,1,1,-1,-1,-1,1, -1,-1,-1,-1,1,1,-1,1, -1,-1,1,1,-1,1,1,-1, 1,1,1,1,1,1,-1,1, 1,1,-1,1,1,-1,-1,1, 1,1,-1,1,-1,-1,-1,1, -1,1,-1,-1,1,-1,-1,1, 1,1,1,1,-1,-1,1,1, -1,-1,1,-1,1,-1,1,1, -1,-1,-1,1,1,-1,-1,-1, -1,1,-1,-1,1,-1,1,1, 1,1,-1,1,-1,1,-1,1, -1,-1,-1,-1,-1,1,-1,1, 1,-1,1,-1,1,1,1,-1, -1,1,-1,-1,-1,1,1,1, -1,-1,-1,-1,-1,-1,-1];
    pilot_screw_rep = repmat(pilot_screw,1,ceil(num_pilots_needed/length(pilot_screw)));
    params.pilots                  = pilot_screw_rep(1:num_pilots_needed).';
    
    params.SC_IND_DATA             = setxor(FILLED_SC_IND,params.SC_IND_PILOTS);     % Data subcarrier indices
    params.N_DATA_IND              = length(params.SC_IND_DATA);             % Number of Data subcarriers
    
    
    % Training Fields
    sts_f = zeros(1,64);
    sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
    sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
    sts_t = ifft(sqrt(13/6).*sts_f, 64);
    sts_t = sts_t(1:16);
    params.sts_t = sts_t;
    
    % LTS for CFO and channel estimation
    % n	N	Preferred Polynomial (1)	Preferred Polynomial (2)
    % 5	31	[5 2 0]	[5 4 3 2 0]
    % 6	63	[6 1 0]	[6 5 2 1 0]
    % 7	127	[7 3 0]	[7 3 2 1 0]
    % 9	511	[9 4 0]	[9 6 4 3 0]
    % 10	1023	[10 3 0]	[10 8 3 2 0]
    % 11	2047	[11 2 0]	[11 8 5 2 0]
    gold_order = log2(params.N_SC);
    switch(params.N_SC)
        case 64
            lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
            params.pilots                  = [1 1 -1 1].';
            params.SC_IND_PILOTS           = [8 22 44 58]; 
            params.SC_IND_DATA             = setxor(FILLED_SC_IND,params.SC_IND_PILOTS); 
        
        case 128
            pol1 = [7 3 0];
            pol2 = [7 3 2 1 0];
            
        case 256
            lts_256 = [-1,-1,-1,1,-1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,1,1,-1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1,1,1,-1,1];
            lts_f(FILLED_SC_IND) = lts_256(1:num_filled_sc);
            
        case 512
            pol1 = [9 4 0];
            pol2 = [9 6 4 3 0];
            
        case 1024
            pol1 = [10 3 0];
            pol2 = [10 8 3 2 0];
            
        case 2048
            pol1 = [11 2 0];
            pol2 = [11 8 5 2 0];
        
        otherwise
            error('Cannot Generate PN Sequence for given N_SC')
    end
    
    if(gold_order~=6 && gold_order~=8)
        init_cond1 = zeros(1,gold_order);
        init_cond2 = zeros(1,gold_order);
        init_cond1(end) = 1;
        init_cond2(end) = 1;
        goldseq = comm.GoldSequence('FirstPolynomial',pol1,...
        'SecondPolynomial',pol2,...
        'FirstInitialConditions',init_cond1,...
        'SecondInitialConditions',init_cond2,...
        'Index',4,'SamplesPerFrame',num_filled_sc);
        gold_pn_seq = goldseq();
        gold_pn_seq = 2*gold_pn_seq - 1;
        
        lts_f = zeros(1,params.N_SC);
        lts_f(FILLED_SC_IND) = gold_pn_seq;
    
    end
        
    
    lts_t = ifft(lts_f, params.N_SC);
    params.lts_f = lts_f;
    params.lts_t = lts_t;


end