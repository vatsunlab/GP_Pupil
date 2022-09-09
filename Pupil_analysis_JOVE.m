function [pupil_data, pupil_dias, PData_mvmt, exclcount,stimorder,uniq_stimSNRs,seqnums] = Pupil_analysis_JOVE(pname,fname)

% re-work of analysis code with hopefully better motion correction

global hFile

% Movement thresholds
mvmt_thresh = 7;  % STDs of diff(MData)
mvmt_time_thresh = 7; % seconds

% choose and load ns2 file and mat file
if nargin==0
    [fname,pname] = uigetfile('D:\DATA\*.ns2','Pick Data File...'); % if single file
end
afname = [fname(1:end-4) '.mat'];
aaa = load(fullfile(pname,afname));
ntrials = numel(aaa.gSM_trial);

[pdata, tdata,NIPStartTime]=ImportStimTrigDatafile(fullfile(pname,fname));
if isempty(pdata)
    exit;
end

pupilidx = find(strcmp('analog 2',{hFile.Entity.Label})); % find the NEV file handle
mvmtidx = find(strcmp('analog 3',{hFile.Entity.Label}));

% determine the stimulus numbers of each SNR
uniq_stimSNRs = unique(aaa.StimCh1.STIMPARAMS(:,8)); 
stimorder = aaa.StimCh1.STIMPARAMS(1,7);
disp([fname ' Order ' num2str(stimorder)]);
noisestimnums = find(aaa.StimCh1.STIMPARAMS(:,8)==0);
stimnums(1) = NaN; seqnums(1) = NaN;
for i = 2:1:numel(uniq_stimSNRs)
    stimnums(i) = find(aaa.StimCh1.STIMPARAMS(:,8)==uniq_stimSNRs(i));
    seqnums(i) = aaa.StimCh1.STIMPARAMS(stimnums(i),8);
end
 
% ---------- DETERMINE MOTION TIMES ----------------%
[~, nsEntityInfo] = ns_GetEntityInfo(hFile, mvmtidx);
recording_length = nsEntityInfo.ItemCount;
 
[~,~,MData] = ...
         ns_GetAnalogData(hFile,mvmtidx,1,recording_length);
     
MData2 = detrend(MData);
std_mvmt = std(MData2);
MData2 = abs(MData2);
MData2(MData2<std_mvmt*mvmt_thresh) = 0;
[~,mvmt_locs] = findpeaks(MData2,1000,'MinPeakDistance',1); % 1000 = sampling rate, 1s between peaks
figure(199); clf; subplot(2,5,[1:3]); 
plot((1:nsEntityInfo.ItemCount)/1000,detrend(MData)); hold on
for i = 1:1:numel(mvmt_locs)
    line([mvmt_locs(i) mvmt_locs(i)],[10 20],'Color','r');
end
set(gca,'YLimMode','manual','YLim',[-20 20]);
ylabel('Motion (a.u.)'); xlabel('Time (s)');
title('Overall session motion');

% extract and average pupil trace around movement times, independent of
% stimuli
PData_mvmt = nan(numel(mvmt_locs),70);
for i = 1:1:numel(mvmt_locs)
   try
    [ns_RES,~,dum] = ...
        ns_GetAnalogData(hFile,pupilidx,mvmt_locs(i)*1000 - 2000,7000); % 2 sec before, to 5 secs after
    catch
        disp('.')
        dum = nan(7000,1);
        ns_RES = 'ns_ERROR';
    end
    if ~strcmpi(ns_RES,'ns_OK')
                disp('.')
                dum = nan(7000,1);
            end
    % downsample, smooth, convert to actual size etc.
    try
        dum = decimate(dum,100); % downsample to 10Hz
    catch
        dum = nan(70,1);
    end
    dum = (dum+5000)*4/10000 + 2; % convert to mm
    
    % subtract 'baseline'
    dum = dum - mean(dum(12:16)); % 500 ms before motion onset
    PData_mvmt(i,1:numel(dum)) = dum; %[0; diff(dum)];
end
figure(199); subplot(2,5,[5 10]); plot((1:size(PData_mvmt,2))/10,nanmean(PData_mvmt,1));
hold on; line([2 2],[-0.1 0.2],'Color','r');
xlabel('Time (s)'); ylabel('{\Delta}PD (mm)');
title('Motion-related pupil change')

% ----------- END DETERMINE MOTION TIMES -------------%

% ----------- Extract full pupil data trace and plot -------------%
[~, nsEntityInfo] = ns_GetEntityInfo(hFile, pupilidx);
recording_length = nsEntityInfo.ItemCount;
 
[~,~,PData] = ...
         ns_GetAnalogData(hFile,pupilidx,1,recording_length);
PData = detrend(PData);
PData = (PData+5000)*4/10000 + 2; % convert to mm

% Identify blinks or large PD changes
PData2 = diff(PData);
[~,blink_locs] = findpeaks(abs(PData2),1000,'MinPeakDistance',0.05,'MinPeakHeight',0.4);

if ~isempty(blink_locs) % blinks detected
    fprintf('Blink correction')
    for kk = 1:1:numel(blink_locs)
        try
        PData(floor(blink_locs(kk)*1e3-100):floor(blink_locs(kk)*1e3-100)+200) = ...
            interp1([-100 100],[PData(floor(blink_locs(kk)*1e3-100)) ...
            PData(floor(blink_locs(kk)*1e3+100))],-100:1:100,'Linear');
        fprintf('.')
        catch
        end
    end
    fprintf('\n');
end
figure(199); subplot(2,5,[6:8]); 
plot((1:nsEntityInfo.ItemCount)/1000,PData); hold on
set(gca,'YLimMode','manual','YLim',[3 5]);
ylabel('Pupil dia. (mm)'); xlabel('Time (s)')
title('Overall session pupil diameter');
% ----------- End plot full pupil trace --------------------------%

pupil_data = {};
pupil_dias = {};
exclcount = zeros(1,numel(uniq_stimSNRs)); pupil_stimtimes = [];

for currstimtype = 1:1:numel(uniq_stimSNRs)
    if currstimtype == 1
        currstimnums = noisestimnums;
    else
        currstimnums = stimnums(currstimtype);
    end
    thisstimtypetrials = [];
    for i = 1:1:numel(currstimnums)
       thisstimtypetrials = [thisstimtypetrials; find(aaa.StimCh1.STIMPARAMS(:,1)==currstimnums(i))];
    end
    thisstimtypetrials = sort(thisstimtypetrials,'ascend');
    thisstimtypetrials(thisstimtypetrials>ntrials-4) = [];
 
    pupil_data{currstimtype} = [];
    pupil_dias{currstimtype} = [];
    
    trialcount = 1; 
    for trialno = 1:1:numel(thisstimtypetrials)
        % find actual stimulus onset - first time trigger goes high
        tdata2=tdata;
        tdata2(tdata2(:,1)<pdata(thisstimtypetrials(trialno),1),:)=[];
        stimstart = tdata2(tdata2(:,2)==32767,:); %all times trigger went high after parallel data       
        stimstart = stimstart(1,1); % first high trigger after parallel data               
        
        dums = stimstart*1e-3 - mvmt_locs; dums(dums<0)=[];
        pupil_stimtimes = [pupil_stimtimes stimstart];
        if any(dums<mvmt_time_thresh)
            Data = zeros(5000,1)./0;
            exclcount(currstimtype) = exclcount(currstimtype)+1;
        else
            
            try
                Data = PData(floor(stimstart-1000):floor(stimstart+3999));
            catch
                Data = nan(5000,1);
            end
            
            %--------- mark stim times on motion and pupil traces ------------%
            figure(199); subplot(2,5,[1:3]);
            if currstimtype==1
                % enabling standard stimulus time markers will slow down
                % program
                %line([stimstart stimstart]*1e-3,[-10 10],'Color',[0.8 0.8 0.8]);
            else
                line([stimstart stimstart]*1e-3,[-10 10],'Color','k');
            end
            figure(199); subplot(2,5,[6:8]);
            if currstimtype==1
                % enabling standard stimulus time markers will slow down
                % program
                %line([stimstart stimstart]*1e-3,[3 5],'Color',[0.8 0.8 0.8]);
            else
                line([stimstart stimstart]*1e-3,[3 5],'Color','k');
            end            
            %-------- end mark stim times on motion and pupil traces --------%            
        end
        
        if isempty(Data) || numel(Data)~=5000
            Data = zeros(5000,1)./0;
        end
        
        % downsample, smooth, convert to actual size etc.
        try
        Data = decimate(Data,100); % downsample to 10Hz
        catch
            Data = zeros(50,1)./0;
        end
        
        try            
            baseline = mean(Data(6:10)); % 0.5 s before stim onset
            Data = (Data - baseline);%/baseline; %change
        catch
            baseline=NaN;
            Data=Data;
        end
        pupil_data{currstimtype}(trialcount,:) = Data;
        pupil_dias{currstimtype} = [pupil_dias{currstimtype}; nanmean(Data(26:33))]; % all pupil dias
        trialcount = trialcount+1;
    end
     disp(['excluded ' num2str(exclcount(currstimtype)) ' trials'])
end



function [parallelData, triggerData, NIPStartTime]=ImportStimTrigDatafile(fname)

global hFile
global fileInfo

[rc, hFile] = ns_OpenFile(fname);
if ~strcmp(rc, 'ns_OK')
    disp('Can''t open this file')
    parallelData = [];
    triggerData = [];
    NIPStartTime = [];
    return;
end
[~,fileInfo] = ns_GetFileInfo(hFile);
nevidx = find(strcmp('NEURALEV',{hFile.FileInfo.FileTypeID})); % find the NEV file handle

parallelID = -1; triggerID = -1;
for eID = 1:1:numel(hFile.Entity(:)),
    [~, entityInfo] = ns_GetEntityInfo(hFile, eID);
    if strcmp(entityInfo.EntityLabel,'Parallel Input'),
        parallelID = eID;
    elseif strcmp(entityInfo.EntityLabel,'SMA 1'),
        triggerID = eID;
    end
end
if triggerID<0 || parallelID<0,
    disp('COULD NOT FIND parallel or trigger inputs');
     parallelData = [];
    triggerData = [];
    NIPStartTime = [];
    return
end

% get stim identifier from parallel in
[~, startIndex] = ns_GetIndexByTime(hFile, parallelID, 0,0);
[~, endIndex] = ns_GetIndexByTime(hFile, parallelID, fileInfo.TimeSpan,0);
nparallelData = 1;
parallelData = zeros(endIndex-startIndex+1,2);
for index=startIndex:endIndex
    [~, ts, eventData, ~] = ns_GetEventData(hFile,parallelID,index);
    eventData = double(eventData);
    parallelData(nparallelData,:) = [ts*1000 eventData];
    nparallelData = nparallelData + 1;
end

% get stim onsets from trigger in
[~, startIndex] = ns_GetIndexByTime(hFile, triggerID, 0,0);
[~, endIndex] = ns_GetIndexByTime(hFile, triggerID, hFile.FileInfo(nevidx).TimeSpan,0);
ntriggerData = 1;
triggerData = zeros(endIndex-startIndex+1,2);
for index=startIndex:endIndex
    [~, ts, eventData, ~] = ns_GetEventData(hFile,triggerID,index);
    eventData = double(eventData);
    triggerData(ntriggerData,:) = [ts*1000 eventData];
    ntriggerData = ntriggerData + 1;
end
NIPStartTime = fileInfo.NIPTime;
