% Generate the noiseless and noisy channel coefficients based on Rayleigh fading
% K :       Number of users
% M :       Number of elements in array antenna
% iters :   Number of trials
% est_SNR : SNR(dB) to estimate the channels

function [H_noisy,H_nonoise]=rayleigh_gen(K,M,iters,est_SNR)
    
    % Noiseless channel
    H_nonoise=zeros(iters,K,M); 
    % Noisy channel
    H_noisy=zeros(iters,K,M);
    
    % Generate noiseless channels
    for iter=1:iters
        H=zeros(K,M);
        for ch=1:K
            H(ch,:)=(normrnd(0,1,[1,M])+1i*normrnd(0,1,[1,M]))/sqrt(2);
        end
        H_nonoise(iter,:,:)=H;
    end
    
    % Noise power
    N_P=-est_SNR;
    N_P=10.^(N_P/10);

    % Add noise
    for iter=1:iters
        noise_vec=zeros(K,M);
        for kk=1:K
        noise_vec(kk,:)=(normrnd(0,sqrt(N_P),[1,M])+1i*normrnd(0,sqrt(N_P),[1,M]))/sqrt(2);
        end
        H_noisy(iter,:,:)=squeeze(H_nonoise(iter,:,:))+noise_vec;
    end
end




