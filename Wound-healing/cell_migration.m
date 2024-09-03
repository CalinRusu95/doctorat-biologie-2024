function [Res_stats,Res_colour,Res_gray,Res_Cells] = cellMigration(dataIn)
    if nargin < 1;     
        help cellMigration; 
        Res_stats = [];
        Res_colour = [];
        Res_gray = [];
        Res_Cells = []; 
        return;  
    end

    if isa(dataIn,'char')
        b = imread(dataIn);    
    else
        b = dataIn;
    end

    if min(b(:)) <  0
        b = b - min(b(:));
    end
    [rows,cols,levs] = size(b);
    f1 = gaussF(rows,cols,1,floor(rows/50),floor(cols/50),1,floor(rows/2),floor(cols/2),1);
    f1 = f1/max(f1(:));
    f2 = gaussF(15,15,1);
    f3 = gaussF(3,3,1);
    
    bf = fftshift(fft2(double(b)));
    if levs>1
        bHP1a = abs(ifft2(((1-f1).*bf(:,:,1))));
        bHP2a = abs(ifft2(((1-f1).*bf(:,:,2))));
        bHP3a = abs(ifft2(((1-f1).*bf(:,:,3))));
        bHP1 = conv2(bHP1a,f2,'same');
        bHP2 = conv2(bHP2a,f2,'same');
        bHP3 = conv2(bHP3a,f2,'same');
        bHP1 = (bHP1-min(bHP1(:)))+(bHP2-min(bHP2(:)))+(bHP3-min(bHP3(:)));
    else
        bHP1a = abs(ifft2(((1-f1).*bf(:,:,1))));
        bHP1 = conv2(bHP1a,f2,'same');
        bHP1 = (bHP1-min(bHP1(:)));
    end

    bHP1 = bHP1(20:end-19,20:end-19);
    [rowsbHP1,colsbHP1] = size(bHP1);
    rowVector = linspace(1,rowsbHP1,2.^floor(log2(rowsbHP1)));
    colVector = linspace(1,colsbHP1,2.^floor(log2(colsbHP1)));
    reduced_bHP1 = interp2(bHP1,colVector,rowVector');
    min_bHP1 = min(reduced_bHP1(:));
    reduced_bHP1 = reduced_bHP1-min_bHP1;
    max_bHP1 = max(reduced_bHP1(:));
    reduced_bHP1 = reduced_bHP1/max_bHP1;
    
    for k = 3:6
        thres2(k-2) = graythresh(reduceu(reduced_bHP1,k));
    end
    thres2 = thres2*max_bHP1+min_bHP1;
    
    thres = min(0.5*thres2);
    
    bCELLS(:,:,1) = (bHP1>(thres));
    angleRange = (0:10:179);
    [TraceCells] = traceTransform_old(bCELLS(1:4:end,1:4:end),angleRange,[-1 1]);
    
    tracesStd(18) = 0;
    for k = 1:18
        ttt = squeeze(TraceCells(:,k,9)); 
        ttt(squeeze(TraceCells(:,k,11))==0) = [];     
        tracesStd(k) = std(ttt);               
    end
    [maxVar,indVar] = max(tracesStd);
    angOrientation = angleRange(indVar);

    sizeOperator = max(10,round(min(rows,cols)/40));
    SE01 = strel('line',sizeOperator,angOrientation);
    SE02 = strel('line',sizeOperator,90+angOrientation);
    
    nhood1                      = getnhood(SE01);
    nhood1                      = (conv2(double(nhood1),gaussF(3,3,1))>0);
    nhood2                      = getnhood(SE02);
    nhood2                      = (conv2(double(nhood2),gaussF(3,3,1))>0);
    SE2                         = strel(nhood2);      
    SE1                         = strel(nhood1);     
    bCELLS2                     = (imclose(bCELLS,SE2));
    bCELLS2                     = (imclose(bCELLS2,SE1));
    bCELLS2a                    = imfill(bCELLS2,'holes');

    bCELLS2a                    = imdilate(bCELLS2a,ones(3));
    if sum(bCELLS2a(:))>(0.95*rows*cols)
        bCELLS3                 = (imopen(bCELLS2,SE2));
    else
        bCELLS3                 = (imopen(bCELLS2a,SE2));
    end
    
    bAreas1                     = bwlabel(bCELLS3);
    bAreas2                     = regionprops(bAreas1,'area');
    bAreas3                     = ismember(bAreas1,find(([bAreas2.Area]/rows/cols)>0.01));
    [bAreas4,numAreas]          = bwlabel(bAreas3);

    bCELLS4 = zerocross(bCELLS3-0.5);
    bCELLS5 = conv2(double(bCELLS4),f3,'same');
    bCELLS6 = bwlabel(bCELLS5>0);
    rp6 = regionprops(bCELLS6,'majoraxislength','perimeter');
    [in1,in2] = sort([rp6.MajorAxisLength]);
    if numAreas==1
        bCELLS7a = bwmorph(bCELLS6==in2(end),'skel',Inf);
    
        bCELLS8a            = padData(bCELLS7a,19,[2 2],1);
        complementArea      = bwlabel(1-bAreas4);
        complementSizes     = regionprops (complementArea,'Area');
        [temp1,temp2]       =max([complementSizes.Area]);
        complementArea2     = (complementArea==temp2);
    
        totArea             = sum(complementArea2(:));
        relArea             = totArea/rows/cols;
        areaCovered         = [totArea;relArea];
        Res_stats.area      = areaCovered;
        commonArea  = padData(complementArea2,19,[2 2],1);
    else
        bCELLS7a            = bwmorph(bCELLS6==in2(end),'skel',Inf);
        bCELLS7b            = bwmorph(bCELLS6==in2(end-1),'skel',Inf);
        bCELLS8a            = padData(bCELLS7a,19,[2 2],1);
        bCELLS8b            = padData(bCELLS7b,19,[2 2],1);
        xyCELLSa            = find(bCELLS8a);
        rowsA               = 1+floor(xyCELLSa/rows);
        colsA               = rem(xyCELLSa-1,rows)+1;
        xyCELLSb            = find(bCELLS8b);
        rowsB               = 1+floor(xyCELLSb/rows);
        colsB               = rem(xyCELLSb-1,rows)+1;
    
        matRA               = repmat(rowsA,[1,size(rowsB,1)]);
        matCA               = repmat(colsA,[1,size(rowsB,1)]);
        matRB               = repmat(rowsB',[size(rowsA,1),1]);
        matCB               = repmat(colsB',[size(rowsA,1),1]);
    
        distBetPoints       = sqrt(((matRA-matRB).^2)+((matCA-matCB).^2  ) );
        [minimumDist1,q1]       = min(distBetPoints); 
        [minimumDist2,q3]       = min(distBetPoints,[],2);
    
        Res_stats.minimumDist   = min(minimumDist1);
        Res_stats.maxDist       = max([max(min(distBetPoints)) max(min(distBetPoints,[],2))   ]);
        Res_stats.avDist        = mean([minimumDist1 minimumDist2']);
    
        Area2A          = bwlabel(1-(conv2(double(bCELLS7a),f3,'same')>0));
        Area2B          = bwlabel(1-(conv2(double(bCELLS7b),f3,'same')>0));
        A1_B1           = sum(sum((Area2A==1)&(Area2B==1)));
        A1_B2           = sum(sum((Area2A==1)&(Area2B==2)));
        A2_B1           = sum(sum((Area2A==2)&(Area2B==1)));
        A2_B2           = sum(sum((Area2A==2)&(Area2B==2)));
        if A1_B1==0
            commonArea  = padData((Area2A==2)&(Area2B==2),19,[2 2],1);
        elseif A1_B2==0
            commonArea  = padData((Area2A==2)&(Area2B==1),19,[2 2],1);
        elseif A2_B1==0
            commonArea  = padData((Area2A==1)&(Area2B==2),19,[2 2],1);
        elseif A2_B2==0
            commonArea  = padData((Area2A==1)&(Area2B==1),19,[2 2],1);
        end
    
        totArea         = sum(commonArea(:));
        relArea         = totArea/rows/cols;
        areaCovered     = [totArea;relArea];
        Res_stats.area  = areaCovered;
    end
    kernelDilation=ones(max(5,ceil(rows/200)));
    
    if levs>1
        if exist('bCELLS8b','var')
            finalBoundaries=uint8(imdilate(bCELLS8a|bCELLS8b,kernelDilation));
        else
            finalBoundaries=uint8(imdilate(bCELLS8a,kernelDilation));
        end
        finalBoundaries = repmat(finalBoundaries,[1 1 3]);
        Res_colour      = 255*finalBoundaries+(1-finalBoundaries).*b;
        Res_colour(:,:,3)=Res_colour(:,:,3).*(1+0.5*uint8(commonArea));
    else
        if exist('bCELLS8b','var')
            finalBoundaries=uint8(imdilate(bCELLS8a|bCELLS8b,kernelDilation));
        else
            finalBoundaries=uint8(imdilate(bCELLS8a,kernelDilation));
        end
        Res_colour      = 255*finalBoundaries+(1-finalBoundaries).*b;
        Res_colour(:,:,2)=Res_colour(:,:,1);
        Res_colour(:,:,3)=0.9*Res_colour(:,:,1).*(1+0.5*uint8(commonArea));
    end
    Res_gray        = double(sum(Res_colour(:,2,:),3)).*(1-0.5*(commonArea));
    
    if nargout==4
        Res_Cells(:,:,1)=bHP1;
        Res_Cells(:,:,2)=bCELLS;
        Res_Cells(:,:,3)=bCELLS2;
        Res_Cells(:,:,4)=bCELLS2a;
        Res_Cells(:,:,5)=bCELLS3;
        Res_Cells(:,:,6)=bCELLS4;
        Res_Cells(:,:,7)=bCELLS5;
        Res_Cells(:,:,8)=bCELLS6;
        Res_Cells(:,:,9)=bCELLS7a;
        if exist('rowsB','var')
            Res_Cells(:,:,10)=bCELLS7b;
        end
    end
    clear a* A* x* t* r* q* n* e* f* i* k* l* m* c* S* T* b* d* P* R1 R2 R3 s*
    