function [lfp_nrem,nrem_psd,rem_psd,lfp_rem,f]=psdd(EEG,sleepstate)
%% REM PSD caculation
n=size(sleepstate.REM,1);rem_psd=cell(1,n);lfp_rem=[];
for i=1:n
    lfprem=EEG(sleepstate.REM(i,1)*1000:sleepstate.REM(i,2)*1000);
    [rempxx,f]=pwelch(lfprem,2048,250,1000,1000);
    rem_psd{i}=rempxx;lfp_rem=[lfp_rem;lfprem];
end
rem_psd=cell2mat(rem_psd);rem_psd=mean(rem_psd,2);
%% NREM PSD caculation
m=size(sleepstate.NREM,1);nrem_psd=cell(1,m);lfp_nrem=[];
for j=1:m
    lfpnrem=EEG(sleepstate.NREM(j,1)*1000:sleepstate.NREM(j,2)*1000);
    [nrempxx,f]=pwelch(lfpnrem,2048,250,1000,1000);
    nrem_psd{j}=nrempxx;lfp_nrem=[lfp_nrem;lfpnrem];
end
nrem_psd=cell2mat(nrem_psd);nrem_psd=mean(nrem_psd,2);