clear;clc;

%讀取檔案
[ppg,ecg,times,number] = textread('PPGECG.txt','%d%d%d%d');

times=times/1000000;


%ppg process -peaks
smoothppg = sgolayfilt(ppg,7,21);
[ppg_pks,ppg_locs,ppg_widths,ppg_proms] = findpeaks(smoothppg,'minpeakheight',450,'MinPeakProminence',20);



ppg_peak_position = size(ppg_locs);
for i = 1 :length(ppg_pks)
ppg_peak_position(i) = times(ppg_locs(i));
end


%ppg process - valleys

vppg = -1*ppg;
smoothvppg = sgolayfilt(vppg,7,21);

[ppg_valleys,ppg_vlocs,ppg_vwidths,ppg_vproms] = findpeaks(smoothvppg,'MinPeakProminence',25,'minpeakheight',-1000);

ppg_valley_position = size(ppg_vlocs);
for i = 1 :length(ppg_valleys)
ppg_valley_position(i) = times(ppg_vlocs(i));
end

%valleys - peaks 比對
ppg_vp_position = size(ppg_valley_position);
ppg_vp = size(ppg_valleys);

for i = 1 :(length(ppg_valleys)-1)
    max = 0;
   for j = 1 :length(ppg_pks) 
      if(ppg_valley_position(i) < ppg_peak_position(j) && ppg_peak_position(j) < ppg_valley_position(i+1))
         if(max<ppg_pks(j))
            max =  ppg_pks(j);
            ppg_vp(i) = ppg_pks(j); 
            ppg_vp_position(i) = ppg_peak_position(j);
         end
      end
      if( ppg_peak_position(j) > ppg_valley_position(i+1))   
          break;
      end
   end
end
%

ppg_record_fre = size(ppg_peak_position-10);
ppg_record_t = size(ppg_peak_position-10);
FREQUENCY = 10;
INDEX = 1;

%算hr
for i = 11:4:(length(ppg_vp_position)-FREQUENCY)
    
                 ppg_record_fre(INDEX) = 600/(ppg_vp_position(i)-ppg_vp_position(i-10));
                 ppg_record_t(INDEX) = ppg_vp_position(i);
                
    INDEX = INDEX + 1 ;
end


%使用常態分佈 算hr

sdppg_record_fre = size(ppg_valley_position-10);
sdppg_record_t = size(ppg_valley_position-10);
sdINDEX = 1;
sd = size(9);

for i = 10:4:(length(ppg_vp_position))
    PIC = 0 ;
    HRtotal = 0;
    for j = 1:9
        sd(j) =  ppg_vp_position (i+1-j) - ppg_vp_position (i-j); 
    
    end
       S = std(sd) ;
       Erd = 3*(S/sqrt(10));
    for k = 1:9
        if(abs(sd(k)-mean(sd))<= Erd ) 
            PIC = PIC + 1 ;
            HRtotal = HRtotal + sd(k) ;
        end
    end
                 sdppg_record_fre(sdINDEX) = (PIC/HRtotal)*60;
                 sdppg_record_t(sdINDEX) = ppg_vp_position(i);
                
    sdINDEX = sdINDEX + 1 ;
end

%fft
%{
fft_ppg = fft(ppg);

figure(10);
plot(fft_ppg);
%}


%ppg部分畫圖
figure(1);
findpeaks(ppg,'minpeakheight',450,'MinPeakProminence',20);
text(ppg_locs+.02,ppg_pks,num2str((1:numel(ppg_pks))'))

figure(11);
plot(times,ppg,'r',times,-smoothvppg,'b')
figure(2);
subplot(2,1,1);plot(ppg_valley_position,-ppg_valleys,'bO',ppg_peak_position,ppg_pks,'rO',ppg_vp_position,ppg_vp,'kO',times,ppg, [14,14],[0,1000],'g-',[44,44],[0,1000],'r-',[63,63],[0,1000],'g-',[92,92],[0,1000],'r-',[107,107],[0,1000],'g-',[137,137],[0,1000],'r-',[162,162],[0,1000],'g-',[201,201],[0,1000],'r-',[222,222],[0,1000],'g-',[247,247],[0,1000],'r-');
title('PPG');
xlabel('S(秒)');
subplot(2,1,2);plot(ppg_record_t,ppg_record_fre, [14,14],[0,350],'g-',[44,44],[0,350],'r-',[63,63],[0,350],'g-',[92,92],[0,350],'r-',[107,107],[0,350],'g-',[137,137],[0,350],'r-',[162,162],[0,350],'g-',[201,201],[0,350],'r-',[222,222],[0,350],'g-',[247,247],[0,350],'r-');
title('PPG');
xlabel('S(秒)');
ylabel('hr');

%ecg process
[ecg_pks,ecg_locs,ecg_widths,ecg_proms] = findpeaks(ecg,'MinPeakProminence',145,'minpeakdistance',15);

ecg_peak_position = size(ecg_locs);
for i = 1 : length(ecg_pks)
ecg_peak_position(i) = times(ecg_locs(i));
end

ecg_record_fre = size(ecg_peak_position-10);
ecg_record_t = size(ecg_peak_position-10);
INDEX = 1;

%算hr

for i = 1:2:(length(ecg_peak_position)-FREQUENCY)
    
                 ecg_record_fre(INDEX) = 600/(ecg_peak_position(i+10)-ecg_peak_position(i));
                 ecg_record_t(INDEX) = ecg_peak_position(i+10);
                
    INDEX = INDEX + 1 ;
end

%ecg部分畫圖
figure(3);
findpeaks(ecg,'MinPeakProminence',145,'minpeakdistance',15);
text(ecg_locs+.02,ecg_pks,num2str((1:numel(ecg_pks))'))

figure(4);
subplot(2,1,1);plot(ecg_peak_position,ecg_pks,'rO',times,ecg, [14,14],[0,1000],'g-',[44,44],[0,1000],'r-',[63,63],[0,1000],'g-',[92,92],[0,1000],'r-',[107,107],[0,1000],'g-',[137,137],[0,1000],'r-',[162,162],[0,1000],'g-',[201,201],[0,1000],'r-',[222,222],[0,1000],'g-',[247,247],[0,1000],'r-');
title('ECG');
xlabel('S(秒)')
subplot(2,1,2);plot(ecg_record_t,ecg_record_fre, [14,14],[0,350],'g-',[44,44],[0,350],'r-',[63,63],[0,350],'g-',[92,92],[0,350],'r-',[107,107],[0,350],'g-',[137,137],[0,350],'r-',[162,162],[0,350],'g-',[201,201],[0,350],'r-',[222,222],[0,350],'g-',[247,247],[0,350],'r-');
title('ECG');
xlabel('S(秒)');
ylabel('hr');


%共同分析
ecg_on_ppg_pks = size(ecg_pks);
for i = 1:length(ecg_peak_position)
    for j = 1:length(times)
        if(ecg_peak_position(i) == times(j))
           ecg_on_ppg_pks(i) = ppg(j);
            break;
        end    
    end                 
end

%共同分析
figure(5);
subplot(2,1,1);plot(ppg_peak_position,ppg_pks,'rO',times,ppg,ecg_peak_position,ecg_pks,'mO',times,ecg, [14,14],[0,1000],'g-',[44,44],[0,1000],'r-',[63,63],[0,1000],'g-',[92,92],[0,1000],'r-',[107,107],[0,1000],'g-',[137,137],[0,1000],'r-',[162,162],[0,1000],'g-',[201,201],[0,1000],'r-',[222,222],[0,1000],'g-',[247,247],[0,1000],'r-');
legend('ppg peak','ppg','ecg peak','ecg');
subplot(2,1,2);plot(sdppg_record_t,sdppg_record_fre,ppg_record_t,ppg_record_fre,ecg_record_t,ecg_record_fre, [14,14],[0,350],'g-',[44,44],[0,350],'r-',[63,63],[0,350],'g-',[92,92],[0,350],'r-',[107,107],[0,350],'g-',[137,137],[0,350],'r-',[162,162],[0,350],'g-',[201,201],[0,350],'r-',[222,222],[0,350],'g-',[247,247],[0,350],'r-');
legend('stdppg hr','ppg hr','ecg hr');

figure(6);plot(ecg_peak_position,ecg_pks+150,'rO',times,ecg+150,ecg_peak_position,ecg_on_ppg_pks-150,'rO',times,ppg-150, [14,14],[0,1000],'g-',[44,44],[0,1000],'r-',[63,63],[0,1000],'g-',[92,92],[0,1000],'r-',[107,107],[0,1000],'g-',[137,137],[0,1000],'r-',[162,162],[0,1000],'g-',[201,201],[0,1000],'r-',[222,222],[0,1000],'g-',[247,247],[0,1000],'r-');
title('ECG & PPG');
xlabel('S(秒)')
grid on;

%ecg 作微分
%smoothECG = sgolayfilt(ECG_data,7,21);
h = 0.01;
difecg = diff(ecg)/h;
GET = 0 ;
THRESHOLD = 3500 ;
CNT = 1 ;

 for i = 1:length(difecg)
        if(difecg(i) > THRESHOLD )
           GET = 1 ;
           difpks(CNT) = difecg(i);
        end  
        if(GET == 1 && difecg(i) < 0  )
            difecg_peak_position(CNT) = times(i);
            difecg_pks(CNT) = difecg(i);
            difecg_on_ecg_pks(CNT) = ecg(i);
           CNT = CNT + 1;
           GET = 0 ;
        end  
 end      

%{
[difecg_pks,difecg_locs,difecg_widths,difecg_proms] = findpeaks(difecg,'MinPeakProminence',7000,'minpeakdistance',30);

difecg_peak_position = size(difecg_locs);
for i = 1 : length(difecg_pks)
difecg_peak_position(i) = times(difecg_locs(i));
end

%ecg 微分圖


figure(7);
findpeaks(difecg,'MinPeakProminence',7000,'minpeakdistance',0);
text(difecg_locs+.02,difecg_pks,num2str((1:numel(difecg_pks))'))




%畫回原ecg上面看結果
difecg_on_ecg_pks = size(difecg_pks);
for i = 1:length(difecg_peak_position)
    for j = 1:length(times)
        if(difecg_peak_position(i) == times(j))
           difecg_on_ecg_pks(i) = ecg(j);
            break;
        end    
    end                 
end
 
%}

figure(8);subplot(2,1,1);plot(times(2:end),difecg,difecg_peak_position,difecg_pks,'rO')
subplot(2,1,2);
plot(ecg_peak_position,ecg_pks+200,'rO',times,ecg+200,difecg_peak_position,difecg_on_ecg_pks-200,'kO',times,ecg-200, [14,14],[0,1000],'g-',[44,44],[0,1000],'r-',[63,63],[0,1000],'g-',[92,92],[0,1000],'r-',[107,107],[0,1000],'g-',[137,137],[0,1000],'r-',[162,162],[0,1000],'g-',[201,201],[0,1000],'r-',[222,222],[0,1000],'g-',[247,247],[0,1000],'r-');
hold on; 
text(difecg_peak_position,difecg_on_ecg_pks-190, num2str((1:numel(difecg_on_ecg_pks))'))
hold off;
