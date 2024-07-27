 function [forces, Pe_step] = LJ_Force(coords,L)
    % Initialize all forces to 0
    forces = zeros(size(coords));

    % Get the number of particles
    nPart = size(coords,2);
    Pe_step = 0;

    % Loop over all particle pairs
    for partA = 1:nPart-1
        for partB = (partA+1):nPart

            % Calculate particle-particle distance
            dr = coords(:,partA) - coords(:,partB);
            % Fix according to periodic boundary conditions
            dr(dr > L/2) = dr(dr > L/2) - L;
            dr(dr < -L/2) = dr(dr < -L/2) + L;
            % Get the distance squared
            dr2 = dot(dr,dr);

            invDr2 = 1.0/dr2; % 1/r^2
            forceFact = invDr2^4 * (invDr2^3 - 0.5);

            % According to Newton's third law, we get action and
            % reaction for the two particles.
            forces(:,partA) = forces(:,partA) + dr*forceFact;
            forces(:,partB) = forces(:,partB) - dr*forceFact;
            if partA ~= partB
                Pe_step = Pe_step + 4*(invDr2^6 - invDr2^3);
            end
        end
    end

    % Multiply all forces by 48
    forces = forces*48;
end