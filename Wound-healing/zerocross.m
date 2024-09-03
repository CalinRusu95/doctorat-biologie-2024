function bor = zerocross(imag);
    
    [lins,cols,levels] = size(imag);
    delta = 0.00001;
    
    if ~isa(imag,'logical')
        imag=sign(imag);
    end
    
    if ((cols==1|lins==1)&levels==1)
        if cols>= lins
            yy = [0 imag(1:cols-1)];
        else
            yy = [0 imag(1:lins-1)']';
        end;
        bor = abs((sign(imag+delta)-sign(yy+delta)));
    elseif (cols==1&lins==1&levels~= 1)
        imag = permute(imag,[2 3 1]);
        yy = [0 imag(1:cols-1)];
        bor = (sign(imag+delta)-sign(yy+delta));
    elseif (lins~= 1&cols~= 1&levels==1)
        diffVer = diff(imag,1,1);zerCols = zeros(1,cols);
        diffHor = diff(imag,1,2);zerRows = zeros(lins,1);
        qq1 = [zerCols;(diffVer)>0];
        qq2 = [(diffVer)<0;zerCols];
        qq3 = [ (diffHor)<0 zerRows ];
        qq4 = [ zerRows (diffHor)>0 ];
        bor = qq1|qq2|qq3|qq4;
    elseif(lins~= 1&cols~= 1&levels~= 1)
        yy5 = [zeros(1,cols,levels);  imag(1:lins-1,1:cols,:)];
        yy6 = [imag(2:lins,1:cols,:);zeros(1,cols,levels)];                  %u|
        yy7 = [ zeros(lins,1,levels) imag(1:lins,1:cols-1,:)];               %-r
        yy8 = [ imag(1:lins,2:cols,:) zeros(lins,1,levels)];                 %l-
        bor5 = fix(delta+(sign(imag)-sign(yy5))/2);
        bor6 = fix(delta+(sign(imag)-sign(yy6))/2);
        bor7 = fix(delta+(sign(imag)-sign(yy7))/2);
        bor8 = fix(delta+(sign(imag)-sign(yy8))/2);
        bor = sign(bor5+bor6+bor7+bor8);
    end
end
