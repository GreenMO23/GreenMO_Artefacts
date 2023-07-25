function [tx_vec_air,loopback,params] = ofdm_tx(params)
%OFDM_TX OFDM Tx samples generation using warplab inspired code.
%Generalized to any number of subcarriers
%  TODO
params.N_DATA_SYMS  = params.N_OFDM_SYMS * length(params.SC_IND_DATA);

sts_t = params.sts_t;
lts_t = params.lts_t;



% Configurable copies of STS, 2.5 Copies of LTS
preamble = [lts_t((params.N_SC/2+1):params.N_SC) repmat(lts_t,1,params.NUM_LTS)];

%% Generate a payload of random integers
if(~params.enable_channel_codes)
    inSig = randi([0 params.MOD_ORDER-1],params.N_DATA_IND,params.N_OFDM_SYMS);
    params.tx_data = inSig(:);
    tx_syms_mat = qammod(inSig,params.MOD_ORDER,'UnitAveragePower',true);
    params.tx_syms_mat = tx_syms_mat;
else
    num_bits_per_sym = log2(params.MOD_ORDER);
    
    if(params.channel_coding_rate==0.5)
        number_of_bits= (params.N_DATA_SYMS * num_bits_per_sym - 2*params.trellis_end_length) * params.channel_coding_rate;
        tx_bits = randi(2, 1, number_of_bits) - 1; 
        % Forward Error Correction 
        tx_bits_coded = double([tx_bits zeros(1,params.trellis_end_length) ]);    % 8 bits padding
        trel = poly2trellis(7, [171 133]);              % Define trellis
        tx_code = convenc(tx_bits_coded,trel);            % convultional encoder
        params.tx_bits = tx_bits;
        params.tx_data = tx_code;
    else
%         puncpat = [1;1;0];
%         K = log2(trel.numInputSymbols); % Number of input streams
%         N = log2(trel.numOutputSymbols); % Number of output streams
%         unpunc_coderate = K/N; % Unpunctured code rate
%         punc_coderate = (K/N)*length(puncpat)/sum(puncpat); % Punctured code rate
%          tx_code = convenc(msg,trellis,puncpat);
        number_of_bits= (params.N_DATA_SYMS * num_bits_per_sym) * params.channel_coding_rate;
        tx_bits = randi(2, number_of_bits,1) - 1;
        convEncoder = comm.ConvolutionalEncoder(poly2trellis(7, [171 133]));
        convEncoder.PuncturePatternSource = 'Property';
        convEncoder.PuncturePattern = [1;1;0;1;1;0];
        tx_code = convEncoder(tx_bits);
        reset(convEncoder)
        params.tx_bits = tx_bits.';
        params.tx_data = tx_code;
    end
    
    
    tx_code_bit_groups = reshape(tx_code,num_bits_per_sym,[]);
    tx_syms = qammod(tx_code_bit_groups,params.MOD_ORDER,'InputType','bit','UnitAveragePower',true);
    mod_bit_size = log2(params.MOD_ORDER);
%   tx_bits = de2bi(tx_data,mod_bit_size);
%     rx_bits_raw = de2bi(op_struct.rx_data(:,pidx),mod_bit_size).';
    tx_syms_mat = reshape(tx_syms, params.N_DATA_IND,params.N_OFDM_SYMS);
    params.tx_syms_mat = tx_syms_mat;
end

% Define the pilot tone values as BPSK symbols
pilots = params.pilots;

% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, params.N_OFDM_SYMS);
params.pilots_mat = pilots_mat;

%% IFFT
% Construct the IFFT input matrix
ifft_in_mat = zeros(params.N_SC, params.N_OFDM_SYMS);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(params.SC_IND_DATA, :)   = tx_syms_mat;
ifft_in_mat(params.SC_IND_PILOTS, :) = pilots_mat;

%Perform the IFFT
tx_payload_mat = ifft(ifft_in_mat, params.N_SC, 1);

% Insert the cyclic prefix
if(params.CP_LEN > 0)
    tx_cp = tx_payload_mat((end-params.CP_LEN+1 : end), :);
    tx_payload_mat = [tx_cp; tx_payload_mat];
end

% Reshape to a vector
tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));
params.tx_payload_vec = tx_payload_vec;

% Construct the full time-domain OFDM waveform
preamble_size = params.num_users*numel(preamble)+(params.num_users)*params.inter_user_zeros;
mimo_preamble_vec = zeros(1,preamble_size);

num_preambs_before = params.user_index; 
one_preamb_size = numel(preamble)+params.inter_user_zeros;
mimo_preamble_vec(num_preambs_before*(one_preamb_size)+1:num_preambs_before*(one_preamb_size)+numel(preamble)) = preamble;

all_preamb_vec = zeros(1,preamble_size);
for usr_idx=1:1: params.num_users
    num_preambs_before = usr_idx-1; 
    one_preamb_size = numel(preamble)+params.inter_user_zeros;
    all_preamb_vec(num_preambs_before*(one_preamb_size)+1:num_preambs_before*(one_preamb_size)+numel(preamble)) = preamble;

end
% mimo_preamble_vec_norm = (0.7 .* mimo_preamble_vec ./ max(abs(mimo_preamble_vec)));
% tx_payload_vec_norm = (0.7 .* tx_payload_vec ./ max(abs(tx_payload_vec)));
% tx_vec = [mimo_preamble_vec_norm tx_payload_vec_norm];
tx_vec = [mimo_preamble_vec tx_payload_vec];
loopback_vec = [all_preamb_vec tx_payload_vec];

% Pad with zeros if required
tx_vec_padded = [tx_vec, zeros(1, params.N_ZERO_PAD)];
loopback_padded = [loopback_vec, zeros(1, params.N_ZERO_PAD)];
tx_vec_air = tx_vec_padded;
loopback_air = loopback_padded;

% Scale the Tx vector to +/- 1
tx_vec_air = (params.TX_SCALE .* tx_vec_air ./ max(abs(tx_vec_air)))*0.8;
loopback_air = (params.TX_SCALE .* loopback_air ./ max(abs(loopback_air)))*0.8;

params.TX_NUM_SAMPS = length(tx_vec_air);
tx_vec_air = tx_vec_air.'; % We deal with columns everywhere else
loopback = loopback_air.'; % We deal with columns everywhere else

% Tx signal
if(params.plot_flag)
    figure(params.cf); clf;

    subplot(2,1,1);
    plot(real(tx_vec_air), 'b');
    axis([0 length(tx_vec_air) -params.TX_SCALE params.TX_SCALE])
    sp1 = gca;
    grid on;
    title('Tx Waveform (I)');

    subplot(2,1,2);
    plot(imag(tx_vec_air), 'r');
    axis([0 length(tx_vec_air) -params.TX_SCALE params.TX_SCALE])
    sp2 = gca;
    grid on;
    title('Tx Waveform (Q)');
    
    linkaxes([sp1 sp2],'x');

    if(params.WRITE_PNG_FILES)
        print(gcf,sprintf('wl_ofdm_plots_%s_txIQ', params.example_mode_string), '-dpng', '-r96', '-painters')
    end
end

end
