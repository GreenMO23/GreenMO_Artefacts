function [channel_struct] = load_channel_struct(meta_params,sim_settings_struct)
    
    %---------------------------------------------------------------------------------
    % Load into local variables for friendlier code 
    num_reflectors = sim_settings_struct.sim_setting.num_reflectors;
    
    
    enable_multipath = sim_settings_struct.enable_multipath;
    if(enable_multipath)
        num_channel_taps = num_reflectors+1; % number of reflectors + direct path
    else
        num_channel_taps = 1; % number of reflectors + direct path
    end
    
    freq = sim_settings_struct.freq;
    wav = 3e8/freq;
    
    multipath_loss = sim_settings_struct.multipath_loss;
    
    reflector_slope_intercepts = sim_settings_struct.sim_setting.reflector_slope_intercepts;
    bts_ant_cords = sim_settings_struct.bts_ant_cords;
    
    user_locs = sim_settings_struct.sim_setting.selected_user_cords;
    [num_users,~] = size(user_locs);
    
    num_ants = sim_settings_struct.num_ants;
    %---------------------------------------------------------------------------------
    % Initialize outputs
    channel_taps_phases = zeros([num_users,num_channel_taps,num_ants]); % in ns
    channel_taps_delays = zeros([num_users,num_channel_taps,num_ants]); % in ns
    channel_taps_mag = zeros([num_users,num_channel_taps,num_ants]);
    channel_vec = zeros([num_users,num_ants]);
    channel_vec_nonoise = zeros([num_users,num_ants]);
    %---------------------------------------------------------------------------------
    % Calculate channel path lengths etc
    for j = 1:1:num_users
      tiled_user_locs = repmat(user_locs(j,:),[num_ants,1]);
      % Norm of order 2, along dimension 2
      tx_rx_dists = vecnorm(tiled_user_locs-bts_ant_cords,2,2).';
      avg_dir_loss = 1/mean(tx_rx_dists);
      if(j==1)
          norm_fac = avg_dir_loss;
      end
      chan_mat = exp(-1j*(2*pi*freq/3e8)*(tx_rx_dists))*avg_dir_loss;
      channel_taps_mag(j,1,:) = avg_dir_loss;
      channel_taps_delays(j,1,:) = mean(tx_rx_dists)/(0.3);
      channel_taps_phases(j,1,:) = exp(-1j*(2*pi*freq/3e8)*(tx_rx_dists));
      
      
      if(enable_multipath)
          for reflector_index=1:1:num_reflectors
            [reflected_distances, path_loss_facs] = get_reflected_dists(tiled_user_locs,bts_ant_cords,reflector_slope_intercepts(reflector_index,:));
            reflected_distances = reflected_distances.';
            avg_refl_loss = 1/(mean(path_loss_facs));
            chan_mat = chan_mat+exp(-1j*(2*pi*freq/3e8)*(reflected_distances))*avg_refl_loss*multipath_loss;
            channel_taps_mag(j,1+reflector_index,:) = avg_refl_loss*multipath_loss;
            channel_taps_delays(j,1+reflector_index,:) = mean(reflected_distances)/(0.3);
            channel_taps_phases(j,1+reflector_index,:) = exp(-1j*(2*pi*freq/3e8)*(reflected_distances));
          end
      end
      % Noise added in reference to user 1
%       analogbeamforming_vec_all_withmp(j,:) = chan_mat;
      
      chan_mat_txpow = chan_mat/norm_fac;
      
      if(j==1)
        noise_pow = 10^((-meta_params.user1_SNR)/10);
        noise_std_dev = sqrt(noise_pow);

      end

      noise_sig = normrnd(0,noise_std_dev/sqrt(2),1,num_ants)+1j*normrnd(0,noise_std_dev/sqrt(2),1,num_ants);
      
      channel_vec_nonoise(j,:) = chan_mat_txpow;
      channel_vec(j,:) = chan_mat_txpow+noise_sig;

  end
  %---------------------------------------------------------------------------------
  % Populate the returning structure
  channel_struct.num_users = num_users;
  channel_struct.max_delay = sim_settings_struct.sim_setting.max_delay;
  channel_struct.num_ants = num_ants;
  channel_struct.num_channel_taps = num_channel_taps;
  channel_struct.taps_mag = channel_taps_mag/norm_fac;
  channel_struct.taps_delays = channel_taps_delays;
  channel_struct.taps_phases =  channel_taps_phases;
  channel_struct.csi_mat = channel_vec;
  channel_struct.csi_mat_nonoise = channel_vec_nonoise;
  channel_struct.samp_freq = meta_params.sampling_freq;
  channel_struct.freq = freq;
  
    
end



function [reflected_dists, path_loss_facs] = get_reflected_dists(tiled_user_locs,dest_pts,reflector_slope_intercept)
    image_pts = zeros(size(dest_pts));
    int_pts = zeros(size(dest_pts));
    if(reflector_slope_intercept(1)~=90)
        m = tan(reflector_slope_intercept(1)*pi/180);
        c = reflector_slope_intercept(2);
        
        x1 = tiled_user_locs(1,1);
        y1 = tiled_user_locs(1,2);
        
        image_pts(:,1) = -2*m*(m*x1-y1+c)/(m^2+1)+x1;
        image_pts(:,2) = 2*(m*x1-y1+c)/(m^2+1)+y1;
        reflected_dists = vecnorm(image_pts-dest_pts,2,2);
        
        m1 = (image_pts(:,2)-dest_pts(:,2))./(image_pts(:,1)-dest_pts(:,1));
        c1 = image_pts(:,2)-m1.*image_pts(:,1);
        
        
        int_pts(:,1) = (c1-c)./(m-m1);
        int_pts(:,2) = m*int_pts(:,1)+c;
        
        path_loss_facs=vecnorm(int_pts-image_pts,2,2).*vecnorm(int_pts-dest_pts,2,2);
        
    else
        %x=c case
        c = reflector_slope_intercept(2);
        
        x1 = tiled_user_locs(1,1);
        y1 = tiled_user_locs(1,2);
        
        image_pts(:,1) = 2*c-x1;
        image_pts(:,2) = y1;
        
        reflected_dists = vecnorm(image_pts-dest_pts,2,2);

        m1 = (image_pts(:,2)-dest_pts(:,2))./(image_pts(:,1)-dest_pts(:,1));
        c1 = image_pts(:,2)-m1.*image_pts(:,1);
        
        int_pts(:,1) = c;
        int_pts(:,2) = m1.*int_pts(:,1)+c1;

        
        path_loss_facs=vecnorm(int_pts-image_pts,2,2).*vecnorm(int_pts-dest_pts,2,2);
    end
end
    
    
    
    
    
