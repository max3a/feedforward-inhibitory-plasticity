function [sleepstate]=stateext(labels)
%% wake state extraction
wake_labels=(labels==2);wake_start=find(diff(wake_labels)>0);wake_stop=find(diff(wake_labels)<0);
if wake_start(1)>wake_stop(1)
   wake_start=[0;wake_start];
end
if wake_start(end)>wake_stop(end)
   wake_stop=[wake_stop;length(labels)];
end
wake=[(wake_start.*2),(wake_stop.*2)];wake_duration=wake(:,2)-wake(:,1);
wake(wake_duration<60,:)=[];
%% REM state extraction
REM_labels=(labels==1);REM_start=find(diff(REM_labels)>0);REM_stop=find(diff(REM_labels)<0);
if REM_start(1)>REM_stop(1)
   REM_start=[0;REM_start];
end
if REM_start(end)>REM_stop(end)
   REM_stop=[REM_stop;length(labels)];
end
REM=[(REM_start.*2),(REM_stop.*2)]; REM_duration=REM(:,2)-REM(:,1);
REM(REM_duration<30,:)=[];
%% NREM state extraction
NREM_labels=(labels==3);NREM_start=find(diff(NREM_labels)>0);NREM_stop=find(diff(NREM_labels)<0);
if NREM_start(1)>NREM_stop(1)
   NREM_start=[0;NREM_start];
end
if NREM_start(end)>NREM_stop(end) 
   NREM_stop=[NREM_stop;length(labels)];
end

NREM=[(NREM_start.*2),(NREM_stop.*2)];NREM_duration=NREM(:,2)-NREM(:,1);
NREM(NREM_duration<60,:)=[];
%% sleep_wake struct
sleepstate=struct('wake',wake,'NREM',NREM,'REM',REM);