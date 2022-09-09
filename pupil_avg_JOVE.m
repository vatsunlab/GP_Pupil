function [datamat,pupil_dia_pct_rl] = pupil_avg_JOVE(animalid, SNR, dB_attn)

% Takes the animalID, SNR, dbAtt as the input

% Generates the final matrix with [animalID, SNR, dbAtt,
% Pupil(1-50)timebinvalues] as the output-1

% Gives the percentage of trials with significant pupil change in that
% session as output-2

%%
% For latin squares design where SNR was not changed in the middle
avg_by_seq = 0;

[fname,pname] = uigetfile('*.ns2','Pick files to average...','Multiselect','On');

if ~iscell(fname)
    fname = {fname};
end
nfiles = numel(fname);

pupil_data = {};
pupil_dias = {};
motion_data = {};
exclcount={};
stimorder = [];
uniq_stimSNRs = {};
% all_stimSNRs = [-24 -18 -12 -6 -3 0 3 6 12 40];
deviant_num = [1 2 3 4 5 6 7 8];
seqnums = {};

for i=1:1:nfiles
    [pupil_data{i}, pupil_dias{i},motion_data{i},exclcount{i}, ...
        stimorder(i),uniq_stimSNRs{i},seqnums{i}] = ...
        Pupil_analysis_JOVE(pname,fname{i});  %Pupil_analysis2(pname,fname{i}); % all individual files%
end


% average same stimulus types
pupil_mean = nan(numel(deviant_num),50); % pupil dia
motion_mean = [];
pupil_std = nan(numel(deviant_num),50);
pupil_dia_means = nan(numel(deviant_num),nfiles);
pupil_dia_bases = nan(1,nfiles);
excl_counts = nan(numel(deviant_num),nfiles);
for i = 1:1:numel(deviant_num)
    dummy1 = []; dummy2 = []; dummy5 = []; dummy6 = [];
    for k = 1:1:nfiles
        
        % average by sequence or by SNR
        if i>1 && avg_by_seq
            index = seqnums{k}(i);
        else
            index = i;
        end
        index=i; %%TEMP ADDED 12/25/2020
        excl_counts(index,k) = exclcount{k}(index);
        
        if i == 1
            motion_mean = [motion_mean; motion_data{k}];
        end
        
        % check if trial complete
        if size(pupil_data{k}{i},2) == 50
            dummy1 = [dummy1; pupil_data{k}{index}];
            
            if i==1
                dummy2 = [dummy2; nanmean(pupil_dias{k}{index})];
                dummy5 = [dummy5; nanstd(pupil_dias{k}{index})./sqrt(numel(pupil_dias{k}{index}))];
            else
                dummy2 = [dummy2; pupil_dias{k}{index}];
                
            end
        else
            dummy1 = [dummy1; zeros(1,50)./0];
            dummy2 = [dummy2; NaN];
            if i==1, dummy5 = NaN; end
        end
    end
    pupil_mean(i,:) = nanmean(dummy1,1);
    pupil_std(i,:) = nanstd(dummy1,[],1)/sqrt(size(dummy1,1)-sum(excl_counts(i,:)));
    pupil_dia_means(i,:) = dummy2;
    
    if i==1, pupil_dia_bases = dummy5; end
        
    switch SNR
        case -24
            SNR = 1;
        case -18
            SNR = 2;
        case -12
            SNR = 3;
        case -6
            SNR = 4;
        case -3
            SNR = 5;
        case 0
            SNR = 6;
        case 3
            SNR = 7;
        case 6
            SNR = 8;
        case 12
            SNR = 9;
        case 40
            SNR = 10;
    end
    
    if i ==1
        SNR = 0;
        temp_data_mat{1,1}= dummy1;
        temp_data_mat{1,2}=[animalid SNR dB_attn];
        temp_data_mat{1,2}= repmat(temp_data_mat{1,2}(1,:),size(dummy1,1),1);
        temp_data_mat{1,3} = [temp_data_mat{1,2},temp_data_mat{1,1}];
    else
        temp_data_mat{2,2}=[animalid SNR dB_attn];
        temp_data_mat{2,1}= vertcat(dummy1,temp_data_mat{2,1});
    end
end
pupil_mean_rl(1,:) = nanmean(temp_data_mat{1,1});
pupil_mean_rl(2,:) = nanmean(temp_data_mat{2,1});
pupil_std_rl(1,:) = nanstd(temp_data_mat{1,1})/sqrt(size(temp_data_mat{1,1},1)-sum(isnan(temp_data_mat{1,1}(:,1))));
pupil_std_rl(2,:) = nanstd(temp_data_mat{2,1})/sqrt(size(temp_data_mat{2,1},1)-sum(isnan(temp_data_mat{2,1}(:,1))));

temp_data_mat{2,2} = repmat(temp_data_mat{2,2}(1,:),size(temp_data_mat{2,1},1),1);
temp_data_mat{2,3} = [temp_data_mat{2,2},temp_data_mat{2,1}];
datamat = vertcat(temp_data_mat{1,3},temp_data_mat{2,3});

disp('Excluded trials...')
disp(sum(excl_counts,2))

% plot overall motion mean and SD
figure;
mm = nanmean(motion_mean,1);

mstd = nanstd(motion_mean,[],1)./sqrt(size(motion_mean,1));
patch([1:1:70 70:-1:1]/10,[mm-mstd fliplr(mm+mstd)],'k','EdgeAlpha',0,'FaceAlpha',0.2);
hold on; plot([1:1:70]/10,nanmean(motion_mean),'Color','k');
plot([2 2],[-0.1 0.25],'Color','r');

title('Motion-linked pupil change');
ylabel('{\Delta}PD (mm)'); xlabel ('Time (s)');

% % % normalize - # of SDs above baseline
pupil_dia_zs = zeros(size(pupil_dia_means));
for i = 1:1:nfiles
    pupil_dia_zs(:,i) = (pupil_dia_means(:,i)-pupil_dia_means(1,i))./pupil_dia_bases(i);
end
pupil_dia_zs(pupil_dia_zs<=2.33) = 0;
pupil_dia_zs(pupil_dia_zs>2.33) = 1; % one-sided test, alpha = 0.01

for i = 1:1:numel(deviant_num)
    
    pupil_dia_pct(i) = nansum(pupil_dia_zs(i,:))/sum(~isnan(pupil_dia_zs(i,:)));
    
end
for i= 1:1:nfiles
    pupil_dia_pct_rl(i) = nansum(pupil_dia_zs(2:end,i))/sum(~isnan(pupil_dia_zs(2:end,i)));
end

pupil_dia_sems = nanstd(pupil_dia_means,[],2)./sqrt(nfiles);

% plot average pupil dia changes

figure; 
%patch([10 20 20 10],[-0.2 -0.2 0.5 0.5],'r','FaceAlpha',0.1,'EdgeAlpha',0); hold on
line([10 10],[-0.2 0.5],'Color','r'); hold on
for i = 1:1:numel(deviant_num)
    if i == 1
        patch([1:1:50 50:-1:1],[pupil_mean(i,:)-pupil_std(i,:) ...
            fliplr(pupil_mean(i,:)+pupil_std(i,:))], ...
            [0 1 0],'EdgeAlpha',0,'FaceAlpha',0.2);
        plot((pupil_mean(i,:)),'Color',[0 1 0]);
    else
        plot((pupil_mean(i,:)),'Color',[10-i 10-i 10-i]./10);
    end
end
line([11 11],[-0.2 0.5],'Color','b');
line([21 21],[-0.2 0.5],'Color','b');
line([35 35],[-0.2 0.5],'Color',[1 0.5 0]);

set(gca,'XLim',[0 50],'YLim',[-0.2 0.5]);
title('Stimulus-linked pupil change')
xticks(0:10:50); xticklabels({'-1','0','1','2','3','4'});
xlabel('Time (s)')
ylabel('{\Delta}PD (mm)')

a=1;