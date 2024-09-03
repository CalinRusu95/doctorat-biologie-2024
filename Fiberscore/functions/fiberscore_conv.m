function [index orientation mean_intensity] = fiberscore_conv(image_name,f,K,L,TC,M,N,T)
   
    disp('Initializing image matrix and parameters...')
    [x,y,D] = size(f);
    index = zeros(x,y);
    
    orientation = ones(x,y)*-100;
    mean_intensity = ones(x,y)*-1;
    
    CMP = zeros(x,y);
    CM = zeros(x,y);
    
    f=double(f)/max(max(double(f))); 
    pad = L+L/2+2;
    P=zeros(x+pad,y+pad);
    
    for i=1:x
        for j=1:y
            P(i+round(pad/2),j+round(pad/2))=f(i,j); 
        end
    end

    
    min_angle = 0; 
    max_angle = (pi/(2*K)*((K))); 
    theta_k = (min_angle:pi/(2*K):max_angle); 
    theta_k_perp = theta_k +pi/2;
    
    mean_fs = cell(size(theta_k));
    nstd_fs = cell(size(theta_k));
    corr_fs = cell(size(theta_k));
    
    mean_fs_perp = cell(size(theta_k_perp));
    nstd_fs_perp = cell(size(theta_k_perp));
    corr_fs_perp = cell(size(theta_k_perp));
    
    I = ((-L/2):(L/2))'; 

    std_dev = 2; 
    mean_profile = 0;
    xx = linspace(-6+mean_profile,6+mean_profile,L+1);
    y_I = exp(-.5*((xx-mean_profile)/std_dev).^2);
    
    for l=1:length(theta_k) 
        
        fprintf('Reading image: %s\n', image_name);
        fprintf('Evaluating the image with theta_k = %.2f degrees\n',theta_k(l)*180/pi);
        fprintf('Angles left to evaluate: %d\n',length(theta_k)-l);
        
        [col row] = pix_displace(theta_k(l),I); 
        [col_perp row_perp] = pix_displace(theta_k_perp(l),I); 
        
        kernel = zeros(L+1); 
        kernel_perp = zeros(L+1);
        kernel_y = zeros(L+1);
        kernel_y_perp = zeros(L+1);
        shift_center = round(length(col)/2); 
        
        for i=1:length(col) 
            kernel(row(i)+shift_center,col(i)*(-1)+shift_center)=1; 
            kernel_perp(row_perp(i)+shift_center,col_perp(i)*(-1)+shift_center)=1; 
            kernel_y(row(i)+shift_center,col(i)*(-1)+shift_center) = y_I(i); 
            kernel_y_perp(row_perp(i)+shift_center,col_perp(i)*(-1)+shift_center) = y_I(i);
        end
        
        term1 = imfilter(f, kernel_y,'conv'); 
        term2 = imfilter(f, kernel,'conv')*sum(y_I)/(L+1);
        term3 =  ((L+1)/L)*(imfilter(f.^2, kernel,'conv')/(L+1)-(imfilter(f, kernel,'conv')/(L+1)).^2);
        term4 = std(y_I);
        std_fs =  stdfilt(f, kernel);    
        mean_fs{l} = imfilter(f, kernel,'conv')*1/(L+1);    
        nstd_fs{l} = std_fs./mean_fs{l};    
        corr_fs{l} = (term1-term2)./(sqrt(term3).*term4*L);
        
        term11 = imfilter(f, kernel_y_perp,'conv'); 
        term22 = imfilter(f, kernel_perp,'conv')*sum(y_I)/(L+1); 
        term33 =  ((L+1)/L)*(imfilter(f.^2, kernel_perp,'conv')/(L+1)- (imfilter(f, kernel_perp,'conv')/(L+1)).^2);
        std_fs_perp =  stdfilt(f, kernel_perp);
        
        mean_fs_perp{l} = imfilter(f, kernel_perp,'conv')*1/(L+1);
        nstd_fs_perp{l} = std_fs_perp./mean_fs_perp{l};
        corr_fs_perp{l} = (term11-term22)./(sqrt(term33).*term4*L);
        
        idxgood = find( corr_fs{l} > TC & nstd_fs{l} > M & nstd_fs{l}./nstd_fs_perp{l} > N & mean_fs_perp{l} > T);
        idxgoodperp = find( corr_fs_perp{l} > TC & nstd_fs_perp{l} > M & nstd_fs_perp{l}./nstd_fs{l} > N & mean_fs{l} > T);
        
        CMP(idxgood) = corr_fs{l}(idxgood); 
        CM(idxgoodperp) = corr_fs_perp{l}(idxgoodperp); 
        
        idxfinal = find(CM > CMP  &  CM > index); 
        idxfinal_perp = find(CMP > CM  &  CMP > index);
        
        index(idxfinal) = CM(idxfinal); 
        orientation(idxfinal) = theta_k(l);
        mean_intensity(idxfinal) = mean_fs{l}(idxfinal);
        
        index(idxfinal_perp) = CMP(idxfinal_perp);
        orientation(idxfinal_perp) = theta_k_perp(l);
        mean_intensity(idxfinal_perp) = mean_fs_perp{l}(idxfinal_perp);
        
    end
end