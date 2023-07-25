function switched_samples = babf2_switcher(meta_params,rx_samples_oversamped)


[~,num_rfc] = size(rx_samples_oversamped);
oversamp_fac = meta_params.oversamp_fac/num_rfc;

clear switched_samples
for rfc_idx=1:1:num_rfc
    curr_samps = resample(rx_samples_oversamped{rfc_idx},1,oversamp_fac);
    undersamp_samples = curr_samps(1+rfc_idx-1:num_rfc:end);
    switched_samples{rfc_idx} = delayseq2(undersamp_samples,(rfc_idx-1)/num_rfc);
end
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins.'));
end
