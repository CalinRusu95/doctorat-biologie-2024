function [shift_i shift_j] = pix_displace(theta,I)

    if (theta <=  45*pi/180) || (theta >=  135*pi/180)
        shift_i = I;
        shift_j = I*(tan(theta)); 
    elseif (theta >  45*pi/180) && (theta <  135*pi/180)
        shift_i = I*(cot(theta)); 
        shift_j = I;
    end
    
    shift_j = round(shift_j);
    shift_i = round(shift_i); 
end