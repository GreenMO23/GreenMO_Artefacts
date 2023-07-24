serial_obj = serial('/dev/ttyUSB12','BaudRate',9600)
fopen(serial_obj);
pause(0.5);
% T.S.             4-3-2-1
% antenna_state =  [[0,1,0,0];...
%                   [1,0,1,0];...
%                   [0,0,1,0];...
%                   [0,1,0,0];...
%                   [1,0,0,0];...
%                   [0,0,0,0];...
%                   [0,0,0,0];...
%                   [0,0,0,0]];
% antenna_decimal_states = antenna_state(:,1)*2^3+antenna_state(:,2)*2^2+antenna_state(:,3)*2+antenna_state(:,4);
antenna_hex_states = '0000000F';
all_antenna_dec =  uint32(hex2dec(antenna_hex_states));
fwrite(serial_obj,all_antenna_dec,'uint32');
pause(0.1)
rx_str_dec = fread(serial_obj, 1, 'uint32');
rx_str_hex = dec2hex(rx_str_dec)
if(all_antenna_dec~=rx_str_dec)
    disp("Antenna state command failed")
end
fclose(serial_obj);
pause(0.5);
clear serial_obj;