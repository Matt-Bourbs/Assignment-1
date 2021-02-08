% Assignment-1
% Matthieu Bourbeau
% 100975211

%% ENHANCEMENTS

    global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons
    
    % Initialization of values and calculation of thermal velocity (v_th).
    particles = 1000;
    region_x = 200e-9;
    region_y = 100e-9;
    step_size = 1e-9;
    timestep = 1000;
    T = 300;
    tau_mn = 0.2e-12;
    v_th = sqrt(2*C.kb*T/C.m_n);
    v_change = step_size/v_th;
    
    % Creates and assigns a random location in the x − y plane.
    % Assigns each particle with a fixed thermal velocity (v_th) and a 
    % random direction.
    % Updates the particle location using Newton's laws of motion at a
    % fixed time interval.
    % Set up a Maxwell-Boltzmann distribution for each velocity component.
    % Check to see if particles are inside boundaries of boxes and if so, move
    % to new location outside boundaries.
    X = rand(2,particles);
    Y = rand(1,particles);
    
    X_position(1,:) = X(1,:)*region_x;
    Y_position(1,:) = Y(1,:)*region_y;
    
    check_X_left = X_position > 0.8e-7;
    check_X_right = X_position < 1.2e-7;
    check_X = check_X_left & check_X_right;
    check_top = Y_position > 0.6e-7;
    check_bottom = Y_position < 0.4e-7;
    box_top = check_top & check_X;
    box_bottom = check_bottom & check_X;
    inside_box = box_top | box_bottom;
    
    while(sum(inside_box) > 0)
        temp_x = rand(1,sum(inside_box));
        temp_y = rand(1,sum(inside_box));
        X_position(inside_box) = temp_x*region_x;
        Y_position(inside_box) = temp_y*region_y;
        
        check_X_left = X_position > 0.8e-7;
        check_X_right = X_position < 1.2e-7;
        check_X = check_X_left & check_X_right;
        check_top = Y_position > 0.6e-7;
        check_bottom = Y_position < 0.4e-7;
        box_top = check_top & check_X;
        box_bottom = check_bottom & check_X;
        inside_box = box_top | box_bottom;
    end
    
    angle(1,:) = X(2,:)*2*pi;

    sigma = sqrt(C.kb*T/C.m_n)/4;
    max_boltz_dist = makedist('Normal',v_th,sigma);
    velocity = random(max_boltz_dist,1,particles);
    figure(1)
    histogram(velocity)
    title('Particle Velocity Histogram')
    X_velocity = v_change*velocity(1,:).*cos(angle(1,:));
    Y_velocity = v_change*velocity(1,:).*sin(angle(1,:));
    
    P_scat = 1-exp(-v_change/tau_mn);
    MFP_vec = zeros(1,particles);
    
    % Add x boundary condition where the particle jumps to the 
    % opposite edge if it reaches the right side it appears at the 
    % left with the same velocity.
    % Add y boundary condition where the particle reflects at the 
    % same angle and retains its velocity.
    % Scatters particles i.e if P_scat is greater than scatter, 
    % the particle scatters.
    % Test to see if particles are inside upper and lower x and y box boundaries and 
    % check if y component must be multiplied by -1.
    for p = 1:timestep
        scatter = rand(1,particles);
        scatter_or_not = scatter < P_scat;
        angle(scatter_or_not) = rand*2*pi;
        
        v = random(max_boltz_dist,1,particles);
        X_velocity(scatter_or_not) = v_change*v(scatter_or_not) ...
            .*cos(angle(scatter_or_not));
        Y_velocity(scatter_or_not) = v_change*v(scatter_or_not) ... 
            .*sin(angle(scatter_or_not));
        
        MFP_vec(~scatter_or_not) = MFP_vec(~scatter_or_not)+step_size;
        MFP_vec(scatter_or_not) = 0;
        
        X_right = X_position + X_velocity > region_x;
        X_position(X_right) = X_position(X_right)+X_velocity(X_right)-region_x;
        
        X_left = X_position + X_velocity < 0;
        X_position(X_left) = X_position(X_left)+X_velocity(X_left)+region_x;
        
        unknown_X = ~(X_left | X_right | box_top | box_bottom);
        X_position(unknown_X) = X_position(unknown_X)+X_velocity(unknown_X);
        
        unknown_Y = (Y_position+Y_velocity > region_y | Y_position+Y_velocity < 0);
        Y_velocity(unknown_Y) = -1*Y_velocity(unknown_Y);
       
        check_X_left = (X_position + X_velocity) > 0.8e-7;
        check_X_right = (X_position + X_velocity) < 1.2e-7;
        check_X = check_X_left & check_X_right;
        check_bottom = (Y_position + Y_velocity) < 0.4e-7;
        box_bottom = check_bottom & check_X;
        X_velocity(box_bottom) = -1*X_velocity(box_bottom);
        
        check_X_left = (X_position + X_velocity) > 0.8e-7 + step_size;
        check_X_right = (X_position + X_velocity) < 1.2e-7 - step_size;
        check_X = check_X_left & check_X_right;
        check_Y_bottom_box_top = Y_position < 0.4e-7 - step_size;
        Y_box_bottom = check_X & check_Y_bottom_box_top;
        Y_velocity(Y_box_bottom) = -1*Y_velocity(Y_box_bottom);
        
        check_X_left = (X_position + X_velocity) > 0.8e-7;
        check_X_right = (X_position + X_velocity) < 1.2e-7;
        check_X = check_X_left & check_X_right;
        check_top = (Y_position + Y_velocity) > 0.6e-7;
        box_top = check_top & check_X;
        X_velocity(box_top) = -1*X_velocity(box_top);
        check_X_left = (X_position + X_velocity) > 0.8e-7 + step_size;
        check_X_right = (X_position + X_velocity) < 1.2e-7 - step_size;
        check_X = check_X_left & check_X_right;
        check_Y_top_box_bottom = Y_position > 0.6e-7 + step_size;
        Y_box_top = check_X & check_Y_top_box_bottom;
        Y_velocity(Y_box_top) = -1*Y_velocity(Y_box_top);
        
        Y_position(1,:) = Y_position(1,:) + Y_velocity(1,:);
        X_store(p,:) = X_position(1,:);
        Y_store(p,:) = Y_position(1,:);
       
        v_avg = sum(velocity)/particles;
        current_temp(1,p) = v_avg^2*C.m_n/(2*C.kb);
        MFP = sum(MFP_vec)/particles;
        TBC_avg = MFP/v_avg;
       
    end
    
    figure(2)
    electron_density = [X_position',Y_position'];
    hist3(electron_density,'CdataMode','auto')
    xlabel('X (m)')
    ylabel('Y (m)')
    title('Electron Density Map')
    colorbar
    view(2)
    
    X_canvas = linspace(0,region_x,10);
    Y_canvas = linspace(0,region_y,10);

    x_bin = discretize(X_position,X_canvas);
    y_bin = discretize(Y_position,Y_canvas);
    temp_bins = zeros(10,10);
    
    for x = 1:10
        for y=1:10
            logicX = x_bin == x;
            logicY = y_bin == y;
            logic = logicX & logicY;
            sum_x = sum(X_velocity(logic))/v_change;
            sum_y = sum(Y_velocity(logic))/v_change;
            v_avg = sqrt((sum_x)^2+(sum_y)^2);
            temp_bins(x,y) = v_avg^2*C.m_n/(2*C.kb);
        end
    end
    
    figure(3)
    surf(temp_bins)
    title('Electron Temperature Map')
    colorbar
    
    box_pic_x = [0.8e-7 0.8e-7  1.2e-7 1.2e-7];
    box_top_y = [1e-7 0.6e-7 0.6e-7 1e-7];
    box_bottom_y = [0 0.4e-7 0.4e-7 0];
    for m = 1:timestep
        figure(4)
        subplot(2,1,1);
        plot(X_store(m,:),Y_store(m,:),'o')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Animated 2-D Particle Trajectories With Box Boundaries')
        hold on
        plot(box_pic_x,box_top_y,'k')
        plot(box_pic_x,box_bottom_y,'k')
        subplot(2,1,2);
        plot(m,current_temp(m),'.k');
        title('Current Temperature Caused')
        ylabel('Temperature (K)')
        xlabel('Time (s)')
        legend(['Current Temperature: ' num2str(current_temp(m))], ...
            ['Avg Time Between Collisions: ' num2str(TBC_avg)], ...
            ['Mean Free Path: ' num2str(MFP)])
        hold on
        xlim([0 timestep])
        ylim([250 350])
        pause(0.00000001)
    end