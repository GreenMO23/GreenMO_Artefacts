function [loopback_noise] = add_noise(loopback,meta_params)
    SNR=meta_params.user1_SNR;
    num_users = meta_params.num_users;
    [num_samps,~] = size(loopback{1});
    for i=1:1:num_users
        first_ltses = loopback{i};
        lts_sensitivity = 10;
        first_ltses = first_ltses(501+lts_sensitivity:500+10.5*64-lts_sensitivity);
        rx_pow = rms(first_ltses)^2;
        noise_pow = rx_pow*(10^((-SNR)/10));
        noise_sig = sqrt(noise_pow/2)*(randn(num_samps,1)+1j*randn(num_samps,1));
        loopback_noise{i} = loopback{i}+noise_sig; 
    end
end

