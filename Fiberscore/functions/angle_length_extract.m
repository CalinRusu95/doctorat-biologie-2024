function [fiber_length orient num max_num_fiber] = angle_length_extract(orient_clean, mean_clean)

    [N M] = size(mean_clean);
    BW = mean_clean ~=0;
    
    [L,num] = bwlabel(BW,8);
    length_label = zeros(N,M);
    
    ii=1;
    max_num_fiber = 0;
    fiber_length = [];
    for i=1:num
        len =sum(sum(L==i));
        id = find(L==i);
        length_label(id) = len;
        if len>=0
            fiber_length(ii) = len;
            ii=ii+1;
        end
        if len> max_num_fiber
            max_num_fiber = len;
        end
    end
    
    orient = orient_clean(BW==1);
    orient = orient(orient~=-100)*180/pi;
end