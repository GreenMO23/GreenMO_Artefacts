function y = apply_channel_to_ofdm(x,fs,h_amp,h_doppler,h_delay)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply a delay-Doppler complex channel to input x
% and output received signal y
% 
% x is a row vector (time domain signal)
% fs is the sampling frequency in Hz
% h_amp is a row vector containing complex gains for each of the paths
% h_doppler is a row vector having the Doppler values in Hertz ( fs/(N*M) )
% h_delay is a row vector having delay values for each path in samples (can be non-inyteger).

% Output:
%       y is of the same size as x after applying channel

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputIsColumn =0;
if(~isrow(x))
    x=x.';
    inputIsColumn =1;
end

Lx = length(x);
nPaths=length(h_amp);
paths = zeros(nPaths,Lx); % variable so that we can have multiple paths with same delay but multiple Dopplers
pathsDelayed = zeros(nPaths,Lx);
for pathIdx = 1:nPaths
    
    paths(pathIdx,:) = h_amp(1,pathIdx)*x.*exp(1i*2*pi*(1/fs)* ...
                                                        h_doppler(pathIdx)*(0:Lx-1)) ;
    pathsDelayed(pathIdx,:) = delayseq(paths(pathIdx,:).', h_delay(pathIdx) ).';
    
%      paths(pathIdx,:) = paths(pathIdx,:) + h_amp(1,pathIdx)*x.*exp(1i*2*pi*(1/fs*4096/Lx)* ...
%                                                          h_doppler(pathIdx)*(0:Lx-1));
% 
%     pathsDelayed(pathIdx,h_delay(pathIdx):(Lx-1)+h_delay(pathIdx)) = paths(pathIdx,:);
end

y = sum(pathsDelayed,1);

if(inputIsColumn)
    y=y.'; % if x was column vector, return a column vector
end
