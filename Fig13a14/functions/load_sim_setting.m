function [sim_params] = load_sim_setting(meta_params, rand_placement)
    
    plot_sim_setting = meta_params.plot_setting;

    sim_params.num_users = meta_params.num_users;    
    sim_params.sim_setting_id = meta_params.setting_id;
    % Load user cords and other settings
    
    sim_params.sim_setting = get_setting_params(meta_params,rand_placement);
    
    sim_params.enable_multipath = true;
    sim_params.freq = 2.4e9;
    
    sim_params.rx_ant_sep_lambda = 0.5;
    sim_params.multipath_loss = 1;
    
    sim_params.wav = 3e8/sim_params.freq;
    sim_params.rx_ant_sep = sim_params.rx_ant_sep_lambda*sim_params.wav;
    
    sim_params.num_ants = meta_params.num_ants;
    sim_params.bts_ant_cords = get_rx_ant_cords(sim_params.sim_setting.bts_cord,meta_params.num_ants, sim_params.rx_ant_sep);
    
    if(plot_sim_setting==1)
        figure(1)
        h=[];
        scat_obj = scatter(sim_params.bts_ant_cords(:,1),sim_params.bts_ant_cords(:,2),100,'s','filled');
        scat_obj.MarkerEdgeColor = "r";
        h(1) = scat_obj;
        hold on
        legend_str= ["BTS"];
        
        color_arr=[[1,0,0];[0,1,0];[0,0,1];[0,1,1];[1,0,1];[1,1,0];[0,0,0];[0.5,1,0.1]];
        for k=1:1:meta_params.num_users
            scat_obj= scatter(sim_params.sim_setting.selected_user_cords(k,1),sim_params.sim_setting.selected_user_cords(k,2),100,'x');
            scat_obj.MarkerEdgeColor = color_arr(k,:);
            h(k+1) = scat_obj;
            legend_str =[legend_str "User Index: "+num2str(k)];
        end

        % Reflector plotting is for setting 1 with vertical and horiz
        % points only
        h(meta_params.num_users+2) = xline(sim_params.sim_setting.grid_dim/2);
        yline(sim_params.sim_setting.reflector_slope_intercepts(1,2));
        xline(-sim_params.sim_setting.grid_dim/2);
        yline(sim_params.sim_setting.reflector_slope_intercepts(4,2));
        
        legend(h,[legend_str 'Reflector']);
        
        h(meta_params.num_users+3) = xline(sim_params.sim_setting.grid_xmin,'--r');
        xline(sim_params.sim_setting.grid_xmax,'--r');
        yline(sim_params.sim_setting.grid_ymin,'--r');
        yline(sim_params.sim_setting.grid_ymax,'--r');
        legend(h,[legend_str 'User location grids']);
        
        ylim([sim_params.sim_setting.plot_ymin,sim_params.sim_setting.plot_ymax]);
        xlim([sim_params.sim_setting.plot_xmin,sim_params.sim_setting.plot_xmax]);
        
        ylim([sim_params.sim_setting.plot_ymin-5,sim_params.sim_setting.plot_ymax]);
        xlim([sim_params.sim_setting.plot_xmin,sim_params.sim_setting.plot_xmax]);
        
    end
end

function [sim_setting_params]= get_setting_params(meta_params,rand_placement)
    setting_id=meta_params.setting_id;
    num_users=meta_params.num_users;
    if(setting_id==1)
       grid_dim = 500;
       sim_setting_params.grid_dim= grid_dim;
       sim_setting_params.max_delay = (2-sqrt(2))*grid_dim/0.3; % in ns
       sim_setting_params.y_offset = grid_dim/10;
       sim_setting_params.bts_cord = [0,grid_dim/5];
       sim_setting_params.grid_ymin = sim_setting_params.bts_cord(2)+sim_setting_params.y_offset;
       sim_setting_params.grid_ymax = grid_dim-5;
       sim_setting_params.grid_xmin = -grid_dim/2+5;
       sim_setting_params.grid_xmax = grid_dim/2-5; % grid for users
       
       sim_setting_params.plot_ymin = -10;
       sim_setting_params.plot_ymax = grid_dim+10;
       sim_setting_params.plot_xmin = -grid_dim/2-10;
       sim_setting_params.plot_xmax = grid_dim/2+10; % plot for the area
       
       sim_setting_params.reflector_slope_intercepts = [[0,grid_dim];[90,grid_dim/2];[90,-grid_dim/2];[0,0]];
        [sim_setting_params.num_reflectors,~] = size(sim_setting_params.reflector_slope_intercepts); 
    end
    % Define more setting ids as you go
    
    if(rand_placement==1)
        geom_spacing = 50;
        sim_setting_params.geom_spacing  = geom_spacing;
        sample_points_vec_x = linspace(sim_setting_params.grid_xmin+2,sim_setting_params.grid_xmax-2,geom_spacing);
        sample_points_vec_y = linspace(sim_setting_params.grid_ymin+2,sim_setting_params.grid_ymax-2,geom_spacing);
        [X,Y] = meshgrid(sample_points_vec_x,sample_points_vec_y);
        sample_points = [X(:) Y(:)];
        num_choices = geom_spacing^2;
        sim_setting_params.selected_user_cords = sample_points(randperm(num_choices,num_users),:);

    % Add code here for debugging specific cases
    elseif(rand_placement==2)
        geom_spacing = 50;
        sim_setting_params.geom_spacing  = geom_spacing;
        sim_setting_params.selected_user_cords  = meta_params.max_users_cords(1:num_users,:);
    elseif(rand_placement==3)
        geom_spacing = 50;
        sim_setting_params.geom_spacing  = geom_spacing;
        sim_setting_params.selected_user_cords  = meta_params.init_users_cords;
    end
    


end


function ant_cords = get_rx_ant_cords(base_station_loc, num_ants, rx_ant_sep, num_rfc, rfc_step_fac)
    if(nargin==3)
        num_rfc=1;
        rfc_step_fac=0;
    end
        
    num_ants_onecol = 1;
    num_ants_onerow = floor(num_ants/num_ants_onecol);
    
    ant_cords = zeros([num_rfc,num_ants,2]);
    
    dist_vec_row = 0:rx_ant_sep:(num_ants_onerow-1)*rx_ant_sep;
    dist_vec_col = 0:rx_ant_sep:(num_ants_onecol-1)*rx_ant_sep;
    rf_chan_sep = num_ants_onerow*rx_ant_sep*rfc_step_fac;
    additional_rfc_incr = rf_chan_sep*[1,0];
    
    for i=1:1:num_rfc
        [X,Y] = meshgrid(dist_vec_row,dist_vec_col);
        cart_prod = [X(:) Y(:)];
        ant_cords(i,:,:) = base_station_loc+(i-1)*additional_rfc_incr+ cart_prod;
    end
    
    if(num_rfc==1)
        ant_cords = squeeze(ant_cords);
    end
end
