function [BW_final] = directional_morph_clean(index,LEN,gdir_thresh, ang_res)
    img = index;
    [x y] = size(img);
     
    min_angle = 0;
    max_angle = min_angle+180;
    
    ang=linspace(min_angle,max_angle,ang_res); 
    D = length(ang);
    union_set= zeros(x,y,D);
    
    for i = 1:D
        se = strel('line',LEN,ang(i));
        lambda = length(se.getneighbors);
        r = lambda;
        if r ==lambda
            IMG = imopen(img,se);
        else
            IMG = ordfilt2(ordfilt2(img,lambda-r+1,se.getnhood)  ,lambda,se.getnhood);
        end
        union_set(:,:,i) = min(img,IMG);
    end
    
    gdir = max(union_set,[],3) - min(union_set,[],3);
    BW_final = gdir>gdir_thresh;
    BW_final = bwmorph(BW_final, 'skel',inf);
    BW_final = bwmorph(BW_final, 'hbreak', inf);
    BW_final = bwmorph(BW_final, 'spur', 5);
end

    
    