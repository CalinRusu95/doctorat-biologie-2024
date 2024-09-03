function [gauss]=gaussF(rowDim,colDim,levDim,rowSigma,colSigma,levSigma,rowMiu,colMiu,levMiu,rho)

    if nargin < 1 
        help gaussF;  
        gauss = []; 
        return; 
    end;
    
    if nargin < 10 
        rho = 0; 
    end;

    if nargin == 1
        [wRow,wCol,wLev] = size(rowDim);
        if wCol == 3      
            levDim = rowDim(3);      
            colDim = rowDim(2);      
            rowDim = rowDim(1);
        elseif wCol == 2  
            colDim = rowDim(2);      
            rowDim = rowDim(1);      
            levDim = 1;
        elseif wCol == 1  
            colDim = 1;      
            levDim = 1;
        end;
    elseif nargin == 2
        levDim = 1;
    end

    filter.x = 1:ceil(rowDim);
    filter.y = 1:ceil(colDim);
    filter.z = 1:ceil(levDim);
    filter.data = zeros(ceil(rowDim),ceil(colDim),ceil(levDim));
    [rr,cc,dd] = meshgrid(filter.x, filter.y ,filter.z);

    if nargin <= 6 
        rowMiu = sum(filter.x)/length(filter.x);
        colMiu = sum(filter.y)/length(filter.y);
        levMiu = sum(filter.z)/length(filter.z);
    end

    sigmVal = 1.1774;
    
    if nargin <= 3    
        rowSigma = (rowMiu-1)/sigmVal;
        colSigma = (colMiu-1)/sigmVal;
        levSigma = (levMiu-1)/sigmVal;
    end;

    rowSigma = max(rowSigma,0.000001);
    colSigma = max(colSigma,0.000001);
    levSigma = max(levSigma,0.000001);
    
    if prod(size(rho)) ~= 1
        if size(rho,1) == 2
            invSigma = inv(rho);
            Srr = invSigma(1,1);
            Scc = invSigma(2,2);
            Src = 2 * invSigma(2,1);
            Srd = 0;
            Scd = 0;
            Sdd = 0;
        else
            invSigma = inv(rho);
            Srr = invSigma(1,1);
            Scc = invSigma(2,2);
            Src = 2*invSigma(2,1);
            Srd = 2*invSigma(1,3);
            Scd = 2*invSigma(2,3);
            Sdd = invSigma(3,3);
        end
        exp_r = (1/rowSigma/rowSigma) * (rr-rowMiu).^2 ;
        exp_c = (1/colSigma/colSigma) * (cc-colMiu).^2 ;
        exp_d = (1/levSigma/levSigma)  *(dd-levMiu).^2;
        exp_rc = (1/rowSigma/colSigma) * (rr-rowMiu).*(cc-colMiu);
        exp_rd = (1/rowSigma/levSigma) * (rr-rowMiu).*(dd-levMiu);
        exp_cd = (1/levSigma/colSigma) * (dd-levMiu).*(cc-colMiu);
        gauss = exp(-(Srr * exp_r + Scc * exp_c  + Sdd * exp_d + Src * exp_rc + Srd * exp_rd + Scd * exp_cd ));
    else 
    
        rho = min(rho,0.999999);
        rho = max(rho,-0.999999);
    
        filter.x2 = (1/(sqrt(2*pi)*rowSigma)) * exp(-((filter.x-rowMiu).^2)/2/rowSigma/rowSigma);
        filter.y2 = (1/(sqrt(2*pi)*colSigma)) * exp(-((filter.y-colMiu).^2)/2/colSigma/colSigma);
        filter.z2 = (1/(sqrt(2*pi)*levSigma)) * exp(-((filter.z-levMiu).^2)/2/levSigma/levSigma);
    
        filter.x2 = filter.x2/sum(filter.x2);
        filter.y2 = filter.y2/sum(filter.y2);
        filter.z2 = filter.z2/sum(filter.z2);
        rhoExponent = (-(rho*(filter.x-rowMiu)' * (filter.y-colMiu))/rowSigma/colSigma);
        filter.rho = (1/sqrt(1-rho^2)) * exp(rhoExponent);
        if (colDim > 1 & rowDim > 1)
            twoDFilter = (filter.x2'*filter.y2).*filter.rho;
            if ceil(levDim) > 1
                for ii = 1:ceil(levDim);
                    threeDFilter(:,:,ii) = twoDFilter.*filter.z2(ii);
                end;
                gauss = threeDFilter;
            else
                gauss = twoDFilter;
            end;
        else 
            if length(filter.x2) > length(filter.y2)
                gauss = filter.x2;
            else
                gauss = filter.y2;
            end
        end;
        gauss(isnan(gauss)) = 0;
    end
end