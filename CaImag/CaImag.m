clear all;
close all;

res = fopen('...\results.txt', 'w+');
specs = fopen('...\spectre.txt', 'w+');

InputPath = '...\Input\...\';

cd(InputPath);
files = dir('**');
files(1:2) = []; 
FilesNum = numel(files); 
cplot = ceil(sqrt(FilesNum));
InSigs = []; 
corrRes = []; 

for i = 1:FilesNum
    
    FileAdrs = strcat(InputPath,files(i).name);
    name = files(i).name;
    opts = detectImportOptions(FileAdrs);
    A = readtable(FileAdrs,opts);
    [M,N] = size(A); 
    
    tstart=360966;
    
    X=A(:,1); 
    spectra=A(:,2:N); 
        
    fracini=0.2;
    fracfin=0.2; 
    pkfrac=0.05; 
    promfrac=0.1; 
    MPD=2000; 
    for j = 1:N-1
        
        spectra = table2array(spectra);
        Yout = smooth(spectra,11,'sgolay'); 
        spectra(:,j)=Yout; 
        intd=gradient(Yout);  
        [c,l] = max(Yout); 
                
        [m,s] = normfit(Yout);
        gaussian = normpdf(Yout,m,s);
        
        X = table2array(X);
        tmax(j)=X(l);
        peakmax(j)=c;
        vmax(j)=max(intd);
        intd_desc = intd(l:M);
        vmin(j)=mean(intd_desc);
        
        first=1;
        last=3;
        
        for k=1:M-1
            if Yout(k+1)>fracini*c && first==1    
                tini(j)=(X(k)+X(k+1))/2;
                kini=k;
                first=2;
            end
            
            if X(k)>tmax(j)
                if Yout(k+1)<fracfin*c && last==3
                    tfin(j)=(X(k)+X(k+1))/2;
                    kfin=k;
                    last=4;
                end;
            end 
        end;
        
        xpeak=X(kini:kfin); 
        peak=Yout(kini:kfin); 
        xpeak2 = typecast(xpeak,'double');
        aria(j)=trapz(xpeak2, peak); 
    
    [pks,locs,widths,proms] = findpeaks(peak,xpeak,'MinPeakHeight',pkfrac*c,'MinPeakDistance',MPD,'MinPeakProminence',promfrac*c);
    figure
	plot(X, spectra)
    title(name)
    
    figure(i)
    subplot(1,2,1)
    plot(X, spectra)
    hold on
    plot(xpeak,peak)
    title(strcat(name,' - signal and peak'))
    subplot(1,2,2)
    plot(xpeak, peak, 'r')
    findpeaks(peak,xpeak,'MinPeakHeight',pkfrac*c,'MinPeakDistance',MPD,'MinPeakProminence',promfrac*c)
    title(strcat(name,' - peak'))
    
    figure(FilesNum+2)
    subplot(cplot,cplot,i)
    plot(X, spectra)
    title(name)
    
    figure(FilesNum+3)
    subplot(cplot,cplot,i)
    plot(X, spectra)
    hold on
    plot(xpeak,peak)
    title(name)
    hold off
    
    Tincrease=tmax-tini; 
    Tdekay=tfin-tmax;
    
    Ttotal=tfin-tini; 
    Latency=tini-tstart; 
    Asymmetry=Tdekay./Tincrease; 
    Increase1Speed=vmax;
    Decrease1Speed=abs(vmin);
    Increase2Speed=0;
    Decrease2Speed=0;
    PeakArea=aria;
    Peak1Amplitude=peakmax; 
    Peak2Amplitude=0;
    PeakNumber = length(pks);
    Type = 1;
    PeaksDistance = 0;
    
    if PeakNumber==2 && pks(1)>pks(2)
        Type = 2;
        PeaksDistance = abs(locs(1)-locs(2));
        Peak1Amplitude = pks(1);
        Peak2Amplitude = pks(2);
        gradis=sort(intd);
        Increase2Speed=gradis(end-1);
        dec_index = find(X==locs(2));
        intd_desc = intd(dec_index:end);
        Decrease2Speed=mean(intd_desc);
    else
        if PeakNumber==2 && pks(1)<pks(2)
            Type = 3;
            PeaksDistance = abs(locs(1)-locs(2));
            Peak1Amplitude = pks(1);
            Peak2Amplitude = pks(2); 
            gradis=sort(intd);
            Increase2Speed=gradis(end-1);
            dec_index = find(X==locs(2));
            intd_desc = intd(dec_index:end);
            Decrease2Speed=mean(intd_desc);
        end
    end
    
    param = [Type, Ttotal, Latency, Asymmetry, Increase1Speed, Increase2Speed, Decrease1Speed, Decrease2Speed, PeakArea, Peak1Amplitude Peak2Amplitude PeaksDistance];
    corrRes(i,1:7)= [Ttotal, Latency, Asymmetry, Increase1Speed, Decrease1Speed, PeakArea, Peak1Amplitude];
    
    fprintf(res, '%9s %3.0f %10.1f %10.1f %15.6f %15.6f %15.6f %15.6f %15.6f %12.4f %11.4f %11.4f %13.0f \n', name, param);
    
    figure(FilesNum+4)
    [H,AX,BigAx,P,PAx] = plotmatrix(corrRes);
    set(gca,'xlim',[0 5])
    ylabel(AX(1,1),'Ttotal','FontSize',10,'Rotation',90)
    ylabel(AX(2,1),'Latency','FontSize',10,'Rotation',90)
    ylabel(AX(3,1),'Asymmetry','FontSize',10,'Rotation',90)
    ylabel(AX(4,1),'Inc. Speed','FontSize',10,'Rotation',90)
    ylabel(AX(5,1),'Dec. Speed','FontSize',10,'Rotation',90)
    ylabel(AX(6,1),'Area','FontSize',10,'Rotation',90)
    ylabel(AX(7,1),'Amplitude','FontSize',10,'Rotation',90)
    
    xlabel(AX(7,1),'Ttotal','FontSize',10,'Rotation',0)
    xlabel(AX(7,2),'Latency','FontSize',10,'Rotation',0)
    xlabel(AX(7,3),'Asymmetry','FontSize',10,'Rotation',0)
    xlabel(AX(7,4),'Inc. Speed','FontSize',10,'Rotation',0)
    xlabel(AX(7,5),'Dec. Speed','FontSize',10,'Rotation',0)
    xlabel(AX(7,6),'Area','FontSize',10,'Rotation',0)
    xlabel(AX(7,7),'Amplitude','FontSize',10,'Rotation',0)

end

fclose(res);
fclose(specs);