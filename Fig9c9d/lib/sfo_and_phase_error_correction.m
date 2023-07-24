function [syms_eq_pc_mat,param] = sfo_and_phase_error_correction(syms_eq_mat,param)
% After channel estimation run this function to apply SFO and phase
% correction

if param.DO_APPLY_SFO_CORRECTION
    % SFO manifests as a frequency-dependent phase whose slope increases
    % over time as the Tx and Rx sample streams drift apart from one
    % another. To correct for this effect, we calculate this phase slope at
    % each OFDM symbol using the pilot tones and use this slope to
    % interpolate a phase correction for each data-bearing subcarrier.
    
    % Extract the pilot tones and "equalize" them by their nominal Tx values
    pilots_f_mat = syms_eq_mat(param.SC_IND_PILOTS, :);
    pilots_f_mat_comp = pilots_f_mat.*param.pilots_mat;
    
    % Calculate the phases of every Rx pilot tone
    pilot_phases = unwrap(angle(fftshift(pilots_f_mat_comp,1)), [], 1);
    
    % Calculate slope of pilot tone phases vs frequency in each OFDM symbol
    pilot_spacing_mat = repmat(mod(diff(fftshift(param.SC_IND_PILOTS)),param.N_SC).', 1, param.N_OFDM_SYMS);
    param.pilot_slope_mat = mean(diff(pilot_phases) ./ pilot_spacing_mat);
    
    % Calculate the SFO correction phases for each OFDM symbol
    pilot_phase_sfo_corr = fftshift((-(param.N_SC/2):(param.N_SC/2-1)).' * param.pilot_slope_mat, 1); 
    pilot_phase_corr = exp(-1i*(pilot_phase_sfo_corr));
    
    % Apply the pilot phase correction per symbol
    syms_eq_mat = syms_eq_mat .* pilot_phase_corr;
else
    % Define an empty SFO correction matrix (used by plotting code below)
    pilot_phase_sfo_corr = zeros(param.N_SC, param.N_OFDM_SYMS);
end


if param.DO_APPLY_PHASE_ERR_CORRECTION
    % Extract the pilots and calculate per-symbol phase error
    pilots_f_mat = syms_eq_mat(param.SC_IND_PILOTS, :);
    pilots_f_mat_comp = pilots_f_mat.*param.pilots_mat;
    pilot_phase_err = angle(mean(pilots_f_mat_comp));
else
    % Define an empty phase correction vector (used by plotting code below)
    pilot_phase_err = zeros(1, param.N_OFDM_SYMS);
end
pilot_phase_err_corr = repmat(pilot_phase_err, param.N_SC, 1);
pilot_phase_corr = exp(-1i*(pilot_phase_err_corr));

param.pilot_phase_err=pilot_phase_err;
% Apply the pilot phase correction per symbol
syms_eq_pc_mat = syms_eq_mat .* pilot_phase_corr;

end