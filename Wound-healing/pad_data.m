function dataPadded=padData (qtdata,numPadPixels,dimsToPad,padWith)
  
    if ~exist('padWith') 
        padWith = 1; 
    end
    
    if (~exist('dimsToPad')) | (isempty(dimsToPad))
        [rows,cols,levs,numFeats] = size(qtdata);
    else
        dimsToPad = [dimsToPad ones(1,4)]; 
        rows = dimsToPad(1);
        cols = dimsToPad(2);
        levs = dimsToPad(3);
        numFeats = dimsToPad(4);
    end
    
    if rows > 1    
        qtdata = [padWith*repmat(qtdata(1,:,:,:),[numPadPixels 1 1]); qtdata;  padWith*repmat(qtdata(end,:,:,:),[numPadPixels 1 1])];
    end

    if cols > 1   
        qtdata = [padWith*repmat(qtdata(:,1,:,:),[1 numPadPixels 1]) qtdata  padWith*repmat(qtdata(:,end,:,:),[1 numPadPixels 1])];
    end
    dataPadded = qtdata;

    if levs > 1
        [rows,cols,levs,numFeats] = size(qtdata);
        qtdata3(:,:,numPadPixels+1:numPadPixels+levs,:)=qtdata;
        qtdata3(:,:,1:numPadPixels,:)=padWith*repmat(qtdata(:,:,1,:),[1 1 numPadPixels]);
        qtdata3(:,:,numPadPixels+levs+1:2*numPadPixels+levs,:)=padWith*repmat(qtdata(:,:,end,:),[1 1 numPadPixels]);
        dataPadded=qtdata3;
    end
end