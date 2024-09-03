function [shift_i, shift_j] = pix_displace(theta, I)
    % Calculate pixel displacements based on the given angle theta
    if isWithinAngleRange(theta, 0, 45) || isWithinAngleRange(theta, 135, 180)
        % For angles close to 0 or 180 degrees, shift primarily in the i-direction
        shift_i = I;
        shift_j = I * tan(theta); 
    else
        % For angles close to 90 degrees, shift primarily in the j-direction
        shift_i = I * cot(theta); 
        shift_j = I;
    end
    
    % Round the shifts to integer values
    shift_i = round(shift_i); 
    shift_j = round(shift_j);
end

% Helper Function

function isWithinRange = isWithinAngleRange(theta, lower, upper)
    % Check if theta is within the specified degree range after converting to radians
    isWithinRange = (theta >= deg2rad(lower)) && (theta <= deg2rad(upper));
end
