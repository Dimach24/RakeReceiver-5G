clc;
clear all; %#ok<CLALL>
close all;
%%
NCellId=17;
caseL   = 'B';
scs     = 30;
pointA  = 4.4;  % GHz
Lmax_   = 8;   % amount of SSB in the HALF-FRAME
mu      = 1;
k_SSB   = 20;
kSSB_bin=int2bit(k_SSB,5,false).';
SFN = 456;
bSFN = int2bit(SFN,10).';
SFN_MSB = bit2int([bSFN(1:6), 0, 0, 0, 0].',10);
SFN_LSB = bit2int(bSFN(7:10).',4);
tran_bandwidth = 60;
toff    =0;
foff    =k_SSB;
samples_offset = 27000;
symbs_received = 60;
MIB     =[...
    0,          ... % just a bit, cos 24 bits required
    bSFN(1:6),   ... % SFN_MSB
    (scs==15||scs==60),     ... % scs15or60
    kSSB_bin(4:-1:1)           ... % kSsbLsb
    1,                      ... % dmrs pos3
    1,1,1,1,                ... % sib.RB=15
    0,1,0,1,    ... % sib.searchSpaceZero=5
    0,          ... % cellBarred=True
    1,          ... % intraFreqReselection=False
    0,          ... % reserved
    ];

%%
samples_gen;
%%
rng(20);

multipath.rays_count=4;
multipath.gains=randn(1,multipath.rays_count-1)*0.2/3+0.4;
multipath.offsets=randi([10,SPS/8],1,multipath.rays_count-1);

%%
[received,copies]=multipath_gen(samples,multipath.rays_count,multipath.offsets,multipath.gains);
%%
r_length=length(received);
match=SsFinder.findPss(received,0,23,SPS);
% match.extra=SsFinder.findPss(copies(1,:),0,23,SPS);

%%

peaks=match.lags(SsFinder.findPeaks(abs(match.corr),0.001));
peaks(peaks<0)=[];
peaks(peaks>peaks(1)+2*SPS)=[];

rake.fingers=length(peaks);
rake.offsets=peaks;
rake.corrs=match.corr(any(match.lags==peaks.',1));
rake.coefs=abs(rake.corrs).^2;
rake.coefs=(rake.coefs/sum(rake.coefs));
rake.copies=zeros(rake.fingers,r_length);

for i=1:rake.fingers
    rake.copies(i,1:r_length-rake.offsets(i)+1)...
    =received(rake.offsets(i):end)*rake.coefs(i);
    % phases=repmat(exp(-1j*2*pi*((1:SPS)-1)/SPS*rake.offsets(i)),1,floor(r_length/SPS));
    % rake.copies(i,:)=ifft(fft(rake.copies(i,:)).*phases);
    % % phase_correction=rake.corrs(i);
    % % phase_correction=phase_correction.'/abs(phase_correction);
    % % angle(phase_correction)/pi*180;
    % % rake.copies(i,:)=rake.copies(i,:)*phase_correction;
end


rake.signal=sum(rake.copies);

rake.match=SsFinder.findPss(rake.signal,0,23,SPS);
rcd=struct();
[rcd.NCellId,rcd.k_SSB,rcd.tindex,rcd.samples]=SsFinder.processSignalByPeakNo(rake.signal,0,23,SPS,1,0.003);
rcd.samples=[rcd.samples, zeros(1,SPS-mod(length(rcd.samples),SPS))];
rcd.rg=OfdmTransceiver.ComplexTime2ResourceGrid(rcd.samples,SPS);
[rcd.pbch,rcd.issb]=ResourceReceiver.getBitstream(rcd.rg,0,rcd.k_SSB,rcd.NCellId,Lmax_);
[rcd.data,rcd.valid_crc]=PbchReceiver.receivePbch(cast(rcd.pbch,"double"),rcd.NCellId,Lmax_);

disp("Rake-receiver status:")
if (rcd.valid_crc)
    disp("data verification success")
else
    disp("data verification failure")
end

[rcd2.NCellId,rcd2.k_SSB,rcd2.tindex,rcd2.samples]=SsFinder.processSignalByPeakNo(received,0,23,SPS,1,0.003);
rcd2.samples=[rcd2.samples, zeros(1,SPS-mod(length(rcd2.samples),SPS))];
rcd2.rg=OfdmTransceiver.ComplexTime2ResourceGrid(rcd2.samples,SPS);
[rcd2.pbch,rcd2.issb]=ResourceReceiver.getBitstream(rcd2.rg,0,rcd2.k_SSB,rcd2.NCellId,Lmax_);
[rcd2.data,rcd2.valid_crc]=PbchReceiver.receivePbch(cast(rcd2.pbch,"double"),rcd2.NCellId,Lmax_);

disp("Regular receiver status:")
if (rcd2.valid_crc)
    disp("data verification success")
else
    disp("data verification failure")
end


%% drawing
figure
[c,l]=xcorr(rake.signal,samples);
c=c/xcorr(rake.signal,rake.signal,0)/xcorr(samples,samples,0);
[~,maxpos]=max(abs(c));
stem(l,abs(c),"Marker","none")
% xlim([l(maxpos)-100,l(maxpos)+100])

%% drawing
figure
subplot(4,1,1)
hold on;
plot(match.lags,abs(match.corr),"Color","b")
% plot(match.lags,abs(match.extra.corr)*multipath.gains(1)/2,"Color","r")
scatter(peaks,abs(match.corr(any(match.lags==peaks.',1))),"Marker","x","SizeData",100,"MarkerEdgeColor","black")
plot(match.lags,zeros(1,length(match.lags))+0.001,"Color","black","LineStyle","--")
xlim([-1000,3e4])
hold off;
title("Корреляция многолучевого сигнала и PSS")
subplot(4,1,2)
plot(rake.match.lags,abs(rake.match.corr),"Color","b")
hold on
plot(rake.match.lags,zeros(1,length(rake.match.lags))+2e-3,"Color","black","LineStyle","--")
title("Корреляция обработанного сигнала и PSS")
xlim([-1000,3e4])
subplot(4,1,4)

rcd2.rg=rcd2.rg./max(abs(rcd2.rg),[],"all");
plt=pcolor(abs(rcd2.rg(1:301,1:end)));
colorbar
plt.EdgeColor='none';
ca=gca();
ca.YDir='normal';
xlim([1,50]);
xlabel('l+1 (номер OFDM символа +1)')
ylabel('k (номер поднесущей)')
title("Восстановленная ресурсная сеть (без RAKE)")


subplot(4,1,3)
rcd.rg=rcd.rg./max(abs(rcd.rg),[],"all");
plt=pcolor(abs(rcd.rg(1:301,1:end)));
colorbar
plt.EdgeColor='none';
ca=gca();
ca.YDir='normal';
xlim([1,50]);
xlabel('l+1 (номер OFDM символа +1)')
ylabel('k (номер поднесущей)')
title("Восстановленная ресурсная сеть (принято RAKE-приемником)")



%% fingers
figure
for i=1:rake.fingers
    [rcdx.NCellId,rcdx.k_SSB,rcdx.tindex,rcdx.samples]=SsFinder.processSignalByPeakNo(received,0,23,SPS,i,0.001);
    rcdx.samples=[rcdx.samples, zeros(1,SPS-mod(length(rcdx.samples),SPS))];
    rcdx.rg=OfdmTransceiver.ComplexTime2ResourceGrid(rcdx.samples,SPS);
    subplot(rake.fingers,1,i)
    rcdx.rg=rcdx.rg./max(abs(rcdx.rg),[],"all");
    plt=pcolor(abs(rcdx.rg(1:301,1:end)));
    colorbar
    plt.EdgeColor='none';
    ca=gca();
    ca.YDir='normal';
    xlim([1,50]);
    xlabel('l+1 (номер OFDM символа +1)')
    ylabel('k (номер поднесущей)')
    title(sprintf("%d finger",i))
end
% %% checking phase errors
% figure
% subplot(3,1,1)
% plt=imagesc(quarter(rcd.rg(1:300,1:end),1e-6));
% colorbar
% % plt.EdgeColor='none';
% ca=gca();
% ca.YDir='normal';
% xlim([1,50]);
% xlabel('l+1 (номер OFDM символа +1)')
% ylabel('k (номер поднесущей)')
% title("четверть (RAKE)")
% subplot(3,1,2)
% plt=imagesc(quarter(rcd2.rg(1:300,1:end),1e-6));
% colorbar
% % plt.EdgeColor='none';
% ca=gca();
% ca.YDir='normal';
% xlim([1,50]);
% xlabel('l+1 (номер OFDM символа +1)')
% ylabel('k (номер поднесущей)')
% title("четверть (без RAKE)")
% 
% subplot(3,1,3)
% plt=imagesc(quarter(rg(1:300,17:end),1e-6));
% 
% colorbar
% % plt.EdgeColor='none';
% ca=gca();
% ca.YDir='normal';
% xlim([1,50]);
% xlabel('l+1 (номер OFDM символа +1)')
% ylabel('k (номер поднесущей)')
% title("четверть (исходный)")
% 
% colormap jet;
%%

fprintf("%.2f%% ошибок (Rake)\n", sum(bits(:,3)~=rcd.pbch.',"all")/length(rcd.pbch)*100);
fprintf("%.2f%% ошибок (без Rake)\n", sum(bits(:,3)~=rcd2.pbch.',"all")/length(rcd2.pbch)*100);