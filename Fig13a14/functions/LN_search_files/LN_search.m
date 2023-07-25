function B_supopt=LN_search(H_noisy,H_nonoise,init_mat)
    
    [K,M] = size(init_mat); 
    % Possible values for a swtich
    All_values=[0 1];
    N=length(All_values);
    
    % nosiy and noiseless channels
    H_correct_channel=H_nonoise;
    
    H_noisy_channel=H_noisy;
    



    % Generate a random initial switch matrix (B)
    B_ini=init_mat;
%         for user_idx=1:K
%             ar = randi([1,N],1,M);
%             B_ini(user_idx,:) = All_values(ar);
%         end

    % Channels
    correct_channel=H_correct_channel;
    noisy_channel=H_noisy_channel;

    % ------ Local search method ------- %

    rate_max=0; % maximum sum-rate
    ok=0; % ok=0 the LN search should be continued
          % ok=1 the LN search should be stopped
    while ok==0
        rate_ini=rate_max; % initial value of sum-rate at each round of the LN searching method
        for B_chan_idx=1:M
            for item=1:N
                for B_user_idx=1:K
                    if B_ini(B_user_idx,B_chan_idx)~=All_values(item)
                        B=B_ini;
                        B(B_user_idx,B_chan_idx)=All_values(item);
                        if sum(B(B_user_idx,:)==0)~=M

                            % received signal
                            received_sig=B*transpose(noisy_channel);

%                             % Noise power
%                             a = 1/sqrt(SNR); 
%                             n_R=(a.*randn(K,M))./sqrt(2);
%                             n_I=(a.*randn(K,M))./sqrt(2);
%                             n=n_R+1i*n_I;
%                             noise_output=B*transpose(n);
%                             sigma2=sum(sum(abs(noise_output).^2));

                            % Obtain SINR
                            SINR=zeros(1,K);
                            for User=1:K
                                Desired_signal=received_sig(:,User);
                                Intrf_signal=received_sig;
                                Intrf_signal(:,User)=[];

                                % Desired power
                                Desired_pow=sum((abs(Desired_signal).^2));

                                % Interference power
                                Intrf_pow=sum(sum((abs(Intrf_signal).^2)));

                                % SINR for each user
                                SINR(1,User)=Desired_pow/(Intrf_pow);
                            end

                            % if B provides a better sum-rate, we
                            % select B as B_opt and update rate_max
                            if sum(log2(1+SINR))>rate_max
                                rate_max=sum(log2(1+SINR));
                                B_opt=B;
                            end
                        end
                    end
                end
            end
        end
        if rate_max==rate_ini
            ok=1;
        else
            B_ini=B_opt;
        end
    end

    % ------ Evaluate based on noiseless channels ------ %
    B=B_opt;

    % Recieved signal
    received_sig=B*transpose(correct_channel);

%     % Noise power
%     a = 1/sqrt(SNR); 
%     n_R=(a.*randn(K,M))./sqrt(2);
%     n_I=(a.*randn(K,M))./sqrt(2);
%     n=n_R+1i*n_I;
%     noise_output=B*transpose(n);
%     sigma2=sum(sum(abs(noise_output).^2));

    % Obtain SINR
    SINR_LN=zeros(1,K);
    for User=1:K
        Desired_signal=received_sig(:,User);
        Intrf_signal=received_sig;
        Intrf_signal(:,User)=[];

        % Desired power
        Desired_pow=sum((abs(Desired_signal).^2));

        % Interference power
        Intrf_pow=sum(sum((abs(Intrf_signal).^2)));

        % SINR for each user
        SINR_LN(1,User)=Desired_pow/(Intrf_pow);
    end

    B_supopt=B_opt;


    
    
end




      
