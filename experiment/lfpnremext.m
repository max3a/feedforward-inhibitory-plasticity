function [lfpnrem,rippleL]=lfpnremext(weak)
%x=bandstop100(EEG);y=bandstop150(x);z=bandstop200(y);m=bandstop250(z);
% rippleL=ripplepassh(m);
%T=(0.001:0.001:length(EEG)/1000)';
%sleepstate=stateext(labels);n=size(sleepstate.NREM,1);lfpnrem=[];%t=[];
%for i=1:n
    %NREM=rippleL(sleepstate.NREM(i,1)*1000:sleepstate.NREM(i,2)*1000,:);
    %time=T(sleepstate.NREM(i,1)*1000:sleepstate.NREM(i,2)*1000,:);
%     lfpnrem=[lfpnrem;NREM];
    %t=[t;time];
%end
% t=(0.001:0.001:length(lfpnrem)/1000);
 %[ripplesL,~,~]=FindRipples([t' lfpnrem]);
 %[~,data] =RippleStats([t' lfpnrem],ripplesL);
 %peakAmp=mean(data.peakAmplitude);duration=mean(data.duration);peakF=mean(data.peakFrequency);incidence=length(data.duration)/length(lfpnrem)*1000;
 %rippleL=struct('peakAmp',peakAmp,'incidence',incidence,'duration',duration,'peakF',peakF);
%end
n=size(weak,2);
for i=1:n
    