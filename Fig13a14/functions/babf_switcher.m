function [switched_samples, filledln_switched_samples,id_switched_samples] = babf_switcher(meta_params,rx_samples_oversamped,channel_est_struct,channel_est_struct_no_noise)

num_rfc=meta_params.num_users;
oversamp_fac = meta_params.oversamp_fac/num_rfc;
num_ants = meta_params.num_ants;
switched_mat = channel_est_struct.babf_switched_mat.';

chan_mat_noisy = channel_est_struct.fully_connected_mat;
chan_mat_noiseless = channel_est_struct_no_noise.fully_connected_mat;

filledln_switched_mat = LN_search(chan_mat_noisy.',chan_mat_noiseless.',channel_est_struct.babf_switched_mat.');

% switched_mat_id = 
clear switched_samples kmeans_switched_samples filledln_switched_samples
for rfc_idx=1:1:num_rfc
    summed_sigs = 0;
    idsummed_sigs = 0;
    filledln_summed_sigs = 0;
    for ant_idx=1:1:num_ants
        curr_samps = resample(rx_samples_oversamped{ant_idx},1,oversamp_fac);
        summed_sigs = summed_sigs + curr_samps*switched_mat(rfc_idx,ant_idx);
        idsummed_sigs = idsummed_sigs + curr_samps*(ant_idx==rfc_idx);
        filledln_summed_sigs = filledln_summed_sigs + curr_samps*filledln_switched_mat(rfc_idx,ant_idx);
    end
    
    undersamp_samples = summed_sigs(1+rfc_idx-1:num_rfc:end);
    idundersamp_samples = idsummed_sigs(1+rfc_idx-1:num_rfc:end);
    filledln_undersamp_samples = filledln_summed_sigs(1+rfc_idx-1:num_rfc:end);
    
    switched_samples{rfc_idx} = delayseq2(undersamp_samples,(rfc_idx-1)/num_rfc);
    id_switched_samples{rfc_idx} = delayseq2(idundersamp_samples,(rfc_idx-1)/num_rfc);
    filledln_switched_samples{rfc_idx} = delayseq2(filledln_undersamp_samples,(rfc_idx-1)/num_rfc);
end
end

function delayed_sig = delayseq2(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ifftshift((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fft(curr_sig);
    delayed_sig = ifft(X .* exp(-1j * frac_del * fbins.'));
end

function filled_ln_mat(chan_mat_noisy,chan_mat_noiseless,init_guess)
    resume=1;
    G1_ini=init_guess; %BABF initialization
    step=0;
    ii=1;
    kk=0;
    S_star_mat=[];
    mu = 0.01;
    divide_rate = 1;
    K = 8;
    M = 64;
    
    while resume==1

        step=step+1;
        jump=0;   % indicator var
        
        % LS_algo_diff_pow ---> local and global search
        [S_star,eval]=LS_algo_diff_pow(iter,K,G1_ini,chan_mat_correct,chan_mat_noisy,1,r);
%         eval1=eval1+eval;
        
        S_star_mat(kk+1,:,:)=S_star;
        
        % obj_func_diff_pow ---> SINR as objective function
        if kk==0 % first time
            min_f=obj_func_diff_pow(iter,K,squeeze(S_star_mat(kk+1,:,:)),chan_mat_correct,chan_mat_noisy);
            % S_final --> optimal switch configuration
            S_final=squeeze(S_star_mat(kk+1,:,:));
        end
        if kk>0
           if obj_func_diff_pow(iter,K,squeeze(S_star_mat(kk+1,:,:)),chan_mat_correct,chan_mat_noisy)>= min_f
               jump=1; %no more processing if bad
           else
               min_f=obj_func_diff_pow(iter,K,squeeze(S_star_mat(kk+1,:,:)),chan_mat_correct,chan_mat_noisy);% update min_f
               S_final=squeeze(S_star_mat(kk+1,:,:));
%                r=r0;
           end
        end

        if jump==0
            ii=1;
            S_star=squeeze(S_star_mat(kk+1,:,:));
            [S_bar,eval]=LS_algo_diff_pow(iter,K,S_star,chan_mat_correct,chan_mat_noisy,2,r);

            eval1=eval1+eval;
            G1_ini = S_bar;
            kk=kk+1;
        end

        if ii>=M*K+1 
            if r>=mu
                r=r/divide_rate;
                ii=1;
                S_star=squeeze(S_star_mat(kk+1,:,:));
                [S_bar,eval]=LS_algo_diff_pow(iter,K,S_star,chan_mat_correct,chan_mat_noisy,2,r);

                eval1=eval1+eval;
                G1_ini = S_bar;
                kk=kk+1;
                jump=0;
            else
                resume=0;
            end
        end    

        if ii<M*K+1 && jump==1

            ii=ii+1;
            kk=kk-1;
            S_star=squeeze(S_star_mat(kk+1,:,:));
            %------
            ii;
            i3=mod(ii-1,M);
            if i3~=0
                i2=((ii-1-i3)/M)+1;
            else
                i2=(ii-1-i3)/M;
                i3=M;
            end
            i2;
            i3;

            S_star_j0=S_star;
            G1_x=zeros(M,1);
            G1_x(i3,1)=1;

            S_star_j0(:,i2)=double(xor(S_star(:,i2),G1_x));
            %------

            [S_bar,eval]=LS_algo_diff_pow(iter,K,S_star_j0,chan_mat_correct,chan_mat_noisy,2,r);
            eval1=eval1+eval;
            G1_ini = S_bar;
            kk=kk+1;
        end
    end
    B_opt_LN_rand=S_final;  % Optimal switch configuration

    % SINR
    Best_SINR=obj_func_diff_pow2(iter,K,B_opt_LN_rand,chan_mat_correct,chan_mat_noisy,SNR_dB);

    SINR_mat(iter,:)=10*log10(Best_SINR);
    waitbar(iter/iters, ff, sprintf('Optimization approach (Progress): %d %%', floor(iter/iters*100)));
end
