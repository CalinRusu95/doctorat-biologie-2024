function [redData] = reduceu(data,numReductions,reducein3D)
  
    if nargin < 1
        help reduceu;   
        redData = [];   
        return; 
    end;
    
    if ~exist('numReductions','var')
        numReductions = 1; 
    end
    
    if ~exist('reducein3D','var')
        reducein3D = 0; 
    end
    
    if ~(isa(data,'double'))
        data = double(data); 
    end
    
    if numReductions == 0
        redData = data;
    else
        if numReductions > 1
            data = reduceu(data,numReductions-1,reducein3D);
        end
    
        [rows,cols,levels] = size(data);
    
        if (levels==1) || (reducein3D == 0)
            redData = convn(data,[1 1;1 1]);
            if levels == 1
                redData  = redData(2:2:end,2:2:end)/4;
            else
                redData  = redData(2:2:end,2:2:end,:)/4;
            end
        else      
            redData = convn(data,ones(2,2,2));
            redData = redData(2:2:end,2:2:end,2:2:end)/8;
    
        end
    end
end