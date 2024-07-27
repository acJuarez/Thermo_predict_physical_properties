function [coords, L] = initCubicGrid(nPart,density)

        % Initialize with zeroes
        coords = zeros(3,nPart);
    
        % Get the cooresponding box size
        L = (nPart/density)^(1.0/3);
    
        % Find the lowest perfect cube greater than or equal to the number of
        % particles
        nCube = 1;
        
        while (nCube^3 < nPart)
            nCube = nCube + 1;
        end
        
        
        % Start positioning - use a 3D index for counting the spots
        index = [0,0,0]';
        
        % Assign particle positions
        for part=1:nPart
            % Set coordinate
            coords(:,part) = (index+[0.5,0.5,0.5]')*(L/nCube);
            
            % Advance the index
            index(1) = index(1) + 1;
            if (index(1) == nCube) 
                index(1) = 0;
                index(2) = index(2) + 1;
                if (index(2) == nCube)
                    index(2) = 0;
                    index(3) = index(3) + 1;
                end
            end
        end
    
    end