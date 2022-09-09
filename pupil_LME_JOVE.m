function pupil_LME_JOVE(datamat)
% takes the final matrix with [animalID, SNR, dbAtt,
% Pupil(1-50)timebinvalues] as the input

M = datamat;
bin_size = 100e-3;
start_tbin = 11; end_tbin = 31;

subjs = unique(M(:,1));
grps = unique(M(:,2));
stepsize = 1/(numel(grps)+1);
colorstmp = [0.85:-stepsize:0]'; colorstmp = (repmat(colorstmp,1,3)); 
colors = [[0 0.447 0.741];colorstmp; ];
figure;
for i = 1:1:numel(grps)
   meandata = nanmean(M(M(:,2)==grps(i),4:end));
   stddata = nanstd(M(M(:,2)==grps(i),4:end))/sqrt(sum(M(:,2)==grps(i)));
   patch([1:numel(meandata) fliplr(1:numel(meandata))]*bin_size, ...
       [meandata-stddata fliplr(meandata+stddata)], colors(i,:), ...
       'FaceAlpha', 0.2, 'EdgeAlpha', 0); hold on
   plot([1:numel(meandata)]*bin_size, meandata, 'Color',colors(i,:),'LineWidth',1); 
end
axis tight
ycoord = get(gca,'YLim');
line([start_tbin start_tbin]*bin_size,[ycoord(1) ycoord(2)],'Color','r','LineStyle','--');
line([end_tbin end_tbin]*bin_size,[ycoord(1) ycoord(2)],'Color','r','LineStyle','--');
xlabel('Time (s)'); ylabel('{\Delta}PD (mm)'); title('Mean {\Delta}PDs')

% compute orthogonal polynomials
 polylen = end_tbin-start_tbin+1; %size(M,2)-5;
tbins = -(polylen-1)/2:(polylen-1)/2; tbins = tbins*bin_size;

%poly1 =  @(x)(x); poly2 =  @(x)(1/2)*(3*x.^2-1); poly3 = @(x) (1/2)*(5*x.^3-3*x);  % Legendre polynomials
poly1 = @(x)(x); poly2 = @(x) (2*x.^2-1); poly3 = @(x) (4*x.^3 - 3*x);% Chebyshev

ot1 = poly1(tbins);
ot2 = poly2(tbins);
ot3 = poly3(tbins);
ot1 = ot1/max(abs(ot1)); ot2 = ot2/max(abs(ot2)); ot3 = ot3/max(abs(ot3));
%figure; plot(tbins,ot1); hold on; plot(tbins,ot2);

% prepare data for GCA
lmedata = nan(polylen*size(M,1),7);

for trials = 1:1:size(M,1)
    lo = (trials-1)*polylen+1;
    hi = trials*polylen;
    lmedata(lo:hi,1) = M(trials,1); %subject
    lmedata(lo:hi,2) = M(trials,2); %group
    lmedata(lo:hi,3) = tbins; % time bin number
    lmedata(lo:hi,4) = ot1; % linear poly
    lmedata(lo:hi,5) = ot2; % quadratic poly
    lmedata(lo:hi,6) = ot3; % cubic poly
    lmedata(lo:hi,7) = M(trials,start_tbin+2:end_tbin+2); % pupil data    
end


lmetab = array2table(lmedata(:,3:7),'VariableNames',{'TBIN','OT1','OT2','OT3','PUPIL'});
lmetab.SUBJ = nominal(lmedata(:,1));
lmetab.GRP = nominal(lmedata(:,2));%,[0:1:numel(grps)]); % specify reference as 0 (standard)
getlabels(lmetab.GRP)


% fit model
lme = fitlme(lmetab,'PUPIL ~ (OT1 + OT2) + GRP + OT1 * GRP + OT2 * GRP + ( 1 | SUBJ )','Verbose',1);
%lme = fitlme(lmetab,'PUPIL ~ (OT1 + OT2 + OT3) + GRP + OT1 * GRP + OT2 * GRP + OT3 * GRP + ( 1 | SUBJ )','Verbose',1);
yfit = fitted(lme);

% plot fitted data averaged across subjects (as subjects have random
% intercept
figure;
for i = 1:1:numel(grps)
   meandata = nanmean(M(M(:,2)==grps(i),3:end));
   stddata = nanstd(M(M(:,2)==grps(i),3:end))/sqrt(sum(M(:,2)==grps(i)));
   errorbar([1:polylen]*bin_size, meandata(start_tbin:end_tbin), ...
       stddata(start_tbin:end_tbin),'Color',colors(i,:), ...
       'LineStyle','none','Marker','o','MarkerFaceColor',colors(i,:), ...
       'MarkerSize',6,'CapSize',0); hold on
end
title('Mean {\Delta}PDs and GCA fits');
xlabel('Time (s)'); ylabel('{\Delta}PD (mm)');

grpav = zeros(numel(grps),polylen);
for grp = 1:1:numel(grps)
    dummy = yfit((lmedata(:,2)==grps(grp)));
    dummy = reshape(dummy,polylen,numel(dummy)/polylen); dummy = dummy';
    grpav(grp,:) = nanmean(dummy);
    plot([1:polylen]*bin_size,grpav(grp,:),'Color','w','LineWidth',3);
    plot([1:polylen]*bin_size,grpav(grp,:),'Color',colors(grp,:),'LineWidth',1);
end
disp(lme)
ncoefs = size(lme.Coefficients.Estimate,1);
contr = eye(ncoefs);
corr_pvals = zeros(ncoefs,1);
for i = 1:1:ncoefs
    corr_pvals(i) = coefTest(lme,contr(i,:),0,'DF','Satterthwaite');
end
disp('Compare original and Satterthwaite-corrected p-vals');
disp([lme.Coefficients.pValue corr_pvals]);

anova(lme,'DF','Satterthwaite')

% maybe a bar of beta values?
aa = reshape(lme.Coefficients.Estimate(4:end),numel(grps)-1,3);
bb = reshape(lme.Coefficients.SE(4:end),numel(grps)-1,3);
cc = reshape(lme.Coefficients.pValue(4:end),numel(grps)-1,3);
dd = reshape(corr_pvals(4:end),numel(grps)-1,3);
cc = dd; % check which p-value to use, dd is Satterthwaite-corrected
cc(cc>0.01) = NaN;cc(cc<=0.01) = 1; 
ngroups = size(aa, 1); nbars = size(aa, 2);

ycoord = max(max(aa))+max(max(bb));

figure; bcols = colormap('lines');bcols(3,:) = [];
for i = 1:1:ngroups
  
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for j = 1:1:3
          
        x = i - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
        bar(x,aa(i,j),0.25,'grouped','EdgeColor','none','FaceColor',bcols(j,:),'FaceAlpha',i/numel(grps)); hold on
        errorbar(x, aa(i,j), bb(i,j), 'k.','CapSize',0);
         text(x,ycoord*cc(i,j),'*','HorizontalAlignment','Center','Color',bcols(j,:));
    end
end
box off
title('GCA weights')
ylabel('{\beta} (a.u.)')
xlabel('SNR (dB)');
xticks(1:1:10);
xticklabels({'-24','-18','-12','-6','-3','0','3','6','12','40'})