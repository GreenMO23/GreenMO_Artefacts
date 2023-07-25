function sum_rate_rand=rand_search(K,M,iters,SNR_dB,H_nonoise)
    
    % Possible values for a swtich
    All_values=[0 1];
    N=length(All_values);
    
    % Noiseless channels
    H_correct_channel=H_nonoise;
    
    % SNR at receiver side
    SNR=10.^(SNR_dB/10);

    % A matrix that saves the rate of each user in different iterations
    rate_rand=zeros(iters,K);

    for iter=1:iters

        % Generate a random initial switch matrix (B)
        B_ini=zeros(K,M);
        for user_idx=1:K
            ar = randi([1,N],1,M);
            B_ini(user_idx,:) = All_values(ar);
        end

        correct_channel=squeeze(H_correct_channel(iter,:,:));

        % ------ Random search method ------- %
        B=B_ini;
        
        % received signal
        received_sig=B*transpose(correct_channel);
        
        % Noise power
        a = 1/sqrt(SNR); 
        n_R=(a.*randn(K,M))./sqrt(2);
        n_I=(a.*randn(K,M))./sqrt(2);
        n=n_R+1i*n_I;
        noise_output=B*transpose(n);
        sigma2=sum(sum(abs(noise_output).^2));

        % Obtain SINR
        SINR_rand=zeros(1,K);
        for User=1:K
            Desired_signal=received_sig(:,User);
            Intrf_signal=received_sig;
            Intrf_signal(:,User)=[];

            % Desired power
            Desired_pow=sum((abs(Desired_signal).^2));

            % Interference power
            Intrf_pow=sum(sum((abs(Intrf_signal).^2)));

            % SINR for each user
            SINR_rand(1,User)=Desired_pow/(sigma2+Intrf_pow);
        end

        % Rate at each user
        rate_rand(iter,:)=log2(1+SINR_rand);
    end
    % Average sum-rate
    sum_rate_rand=sum(mean(rate_rand));

end




      
