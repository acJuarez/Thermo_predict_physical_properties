clear all    % clear all variables from the Worespace
close all    % close all open figures
clc          % clear Command Window

% ===================
%     Initialize
% ===================

% Set configuration parameters
nPart = 5^3;        % Number of particles - ideally a cubic number  
density = 0.05;     % Density of particles
nSteps = 2500;     % Total simulation time (in integration steps)
dt = 0.005;        % Integration time
Temp = [1.8 1.9 2, 2.1 2.2]; % Simulation temperature
TotalE = zeros(5,1);
eTemp = zeros(5,1);
kb = 1.380649*10^(-23);
for j = 1:5
    
    % Set Andersen thermostat parameters
    thermostatFlag = 0; % set to 1 to turn it on
    nu = 1;          % frequency of collisions with heatbath
    velSTD = Temp(j)^0.5;  % the standard deviation for velocities drawn from the Boltzmann distribution at some temperature
    
    % Set initial configuration
    [position_t, L] = initCubicGrid(nPart, density); 
    
    % Initialize arrays to store information from the simulation
    position_all = zeros(3, nPart, nSteps+1);
    velocity_all = zeros(3, nPart, nSteps+1);
    acceleration_all = zeros(3, nPart, nSteps+1);
    displacement_all = zeros(3, nPart, nSteps+1);
    MSD = zeros(1, nSteps);
    Ke_t = zeros(nSteps +1,1);
    Pe_t = zeros(nSteps +1,1);
    
    %% Assignment 3.a: initialize velocities
    % set RNG seed
    seed = 12345;
    rng ( seed );
    
    % generate uniform random numbers in the range [0, 1] for each velocity component
    velocity_t = rand(size(position_t));
    
    % net translation = 0
    velocity_t = velocity_t - mean(velocity_t, 2);
    
    % scale to temperature
    totV2 = sum(velocity_t.^2, 'all')/nPart;  % Mean-squared velocity
    velScale = sqrt(3*Temp(j)/totV2);   % Velocity scaling factor
    velocity_t = velocity_t*velScale; 
    
    %% Assignment 3.f: calculate initial forces
    
    % Calculate initial accelerations & initial PE
    [acceleration_t, Pe_step] = LJ_Force(position_t, L);
    Pe_t(1) = Pe_step;
    %acceleration_t = 0;
    
    %%
    % save initial information
    position_all(:,:,1) = position_t;
    velocity_all(:,:,1) = velocity_t;
    acceleration_all(:,:,1) = acceleration_t;
        for i = 1 : nPart
            vel_mag = sqrt(velocity_t(1,i)^2 + velocity_t(2, i)^2 + velocity_t(3, i)^2);
            Ke_t(1) = Ke_t(1) + 0.5*(vel_mag).^2;
        end
    
    %%
    % ===================
    % Molecular Dynamics Loop
    % ===================
    for step = 1:nSteps
        % showAnimation(position_t, L)  % this function is provided
        % Comment out the above line to speed up performance by not visualizing
        % the particles
        
        % First Verlet update: position
        dr = dt*velocity_t + 0.5*dt^2*acceleration_t;
        displacement_all(:,:,step+1) = dr+displacement_all(:,:,step); % store displacement from original position
        MSD(step) =  1/nPart * sum( sum(displacement_all(:,:,step+1).^2) );
        position_t = position_t + dr;
    
        % correct positions for PBC
        position_t(position_t > L) = position_t(position_t > L) - L;
        position_t(position_t < 0) = position_t(position_t < 0) + L;
        
        % store particle positions
        position_all(:,:,step+1) = position_t;
        
        % First velocity update with the current acceleration
        velocity_t = velocity_t + 0.5*acceleration_t*dt;
        
        %% Assignment 3.e: calculate accelerations
        [acceleration_t, Pe_step] = LJ_Force(position_t, L);
        
        %%
        % Store new acceleration
        acceleration_all(:,:,step+1) = acceleration_t;
        Pe_t(step+1) = Pe_step;
        
        % Second velocity update with the new acceleration
        velocity_t = velocity_t + 0.5*acceleration_t*dt;
        
        % store particle velocities & kinetic energy 
        velocity_all(:,:,step+1) = velocity_t;
        for i = 1 : nPart
            vel_mag = sqrt(velocity_t(1,i)^2 + velocity_t(2, i)^2 + velocity_t(3, i)^2);
            Ke_t(step+1) = Ke_t(step+1) + 0.5*(vel_mag).^2;
        end
        
    
        
        % Implement the Andersen thermostat
        if thermostatFlag
            for part =1:nPart
                % Test for collisions with the Andersen heat bath
                if (rand < nu*dt)
                    % If the particle collided with the bath, draw a new velocity
                    % out of a normal distribution
                    velocity_t(:,part) = normrnd(0,velSTD,3,1);
                end
            end
        end
            
        if mod(step,100) == 0
            fprintf('Step %i of %i \n', step, nSteps);
        end
    end
    
    figure
    title("Plot of particle trajectories")
    for particle = 1:10
        scatter3(displacement_all(1,particle,:), displacement_all(2,particle,:), displacement_all(3,particle,:))
        hold on
    end
    hold off
    
    %calculate equilibrium temp
   
    eTemp(j) =(mean(Ke_t))*(2/3)*(1/nPart);
    
    %Plot for kinetic, potential, and total energy
    TotalE(j) = mean(Ke_t + Pe_t);
    t = linspace(0, nSteps, nSteps+1);
    figure
    title('Plot of Kinetic, Potential, & Total Energy')
    plot(t, Ke_t);
    hold on
    plot(t, Pe_t);
    plot(t, Pe_t + Ke_t);
    legend('Kinetic Energy', 'Potential Energy', 'Total Energy');
    hold off



%% 4.D. From looking at graph, long-time domain begins at about nstep 400
% t = 1:nSteps;
% MSD = t.^2;
% %set long domain boundary at about 400 for tests, 1000 for full run
% linRegStart = 1000;
% linReg = fitlm( table(t(linRegStart:end)', MSD(linRegStart:end)'),'linear');
% 
% figure();
% hold on;
% plot(t, MSD, 'k', 'linewidth', 2)
% plot(linReg)
% xlabel("nSteps");
% ylabel("Mean-square displacement");
% % legend("Total energy", "location", "east");
% title('Mean square displacement against time')

%% 4.e linear regression for five eTemp
figure();
hold on;
linReg_5 = fitlm( table(eTemp, TotalE),'linear')
plot(linReg_5);
title('Total energy against temperature');
xlabel('Dimensionless equilibrium temperature, T_eq');
ylabel('Dimensionless total internal energy, U_tot')
end 