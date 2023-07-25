clc
clear 
close all

% Number of users
K=4;
% Number of elements in the array antenna
M=32;
% SNR of the channel estimation process
SNR_est=20;
% SNR(dB) at the receiver
SNR=20;

% Numer of trials
iters=10;

% noise_cond=false : Noiseless channel are used in our algorithm  
% noise_cond=true  : Noisy channel are used in our algorithm   
noise_cond=true; 

% Generate the channels
[H_noisy,H_nonoise]=rayleigh_gen(K,M,iters,SNR_est);

% sum-rate with local neighbor search
sum_rate_LN=LN_search(K,M,iters,SNR,H_noisy,H_nonoise,noise_cond)

% sum-rate with random search
sum_rate_rand=rand_search(K,M,iters,SNR,H_nonoise)

