function [td_channel] = plot_channel_estimates(rx_H_est,param,os_fac)
% plot rx_H_est (db and angle) in time and frequency domain
% Input
%     rx_H_est: freq-time channel (not fftshifted)
%     param: not needed: only for backward-compatibility

% rx_H_est = rx_H_est_vec(pack_ind,:);
if(nargin<3)
    os_fac = 1;
end

if(nargin<2)
    param.N_SC = length(rx_H_est);
end

subc_ind = (1:param.N_SC) - param.N_SC/2 - 1;
td_channel = ifft(rx_H_est,os_fac*param.N_SC);
fd_channel = rx_H_est;

sp1=subplot(221);
plot(subc_ind,db(fftshift(fd_channel)))
hold on; grid on;

set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([-32 31])
% ylim([4 12])
xlabel('Subcarrier Index','FontSize',12)
ylabel('CSI (dB)','FontSize',12)
title('CSI','FontSize',12)


sp3=subplot(223);
plot(subc_ind,unwrap(angle(fftshift(fd_channel))))
hold on; grid on;

set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([-10 -6])
xlabel('Subcarrier Index','FontSize',12)
ylabel('Phase','FontSize',12)

sp2 = subplot(222);
plot(abs(td_channel));
hold on; grid on;

set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('Sample Idx','FontSize',12)
ylabel('CIR (dB)','FontSize',12)
% xlim([1 64])
title('CIR','FontSize',12)

sp4 = subplot(224);
plot(angle(td_channel));
hold on; grid on;

set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('Sample Idx','FontSize',12)
ylabel('Phase','FontSize',12)
% xlim([1 64])

linkaxes([sp1 sp3],'x');
linkaxes([sp2 sp4],'x');


