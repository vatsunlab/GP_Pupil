function pupil_threshold_estimate_JOVE(pupil_dia_pct_mat)

% Input-1
% pupil_dia_pct_mat: Put all the session-wise percentage of trials with significant 
% pupil changes into each cell of a cell array where the cells are arranged from
% lower to higher SNR

% code outputs psychometric curve and estimates threshold 

%%
M1 = pupil_dia_pct_mat; 

all_stimSNRs = [-24, -18, -12, -6, -3, 0, 3, 6, 12, 40];
for i = 1:1:size(M1,2)
    dims(1,i) = numel(M1{i});
end
dim1 = max(dims);
dim2 = size(M1,2);
pupil_dia_pct = NaN(dim1,dim2);


for i = 1:1:size(M1,2)
    for j = 1:1:numel(M1{i})
        pupil_dia_pct(j,i) = M1{i}(j);
    end
end

for i=1:1:size(pupil_dia_pct,2)
    pupil_pct_std(1,i) = nanstd(pupil_dia_pct(:,i));
    pupil_pct_sem(1,i) = pupil_pct_std(1,i)./sum(~isnan(pupil_dia_pct(:,i)));
end  


statset('Robust','on');
data1 = pupil_dia_pct;
SNRs = repmat(all_stimSNRs,size(data1,1),1);
SNR = SNRs(:);
fitsnrs = -40:0.1:40;

% fit to percent trials on which significant change was observed
y2 = pupil_dia_pct(:);
sigmf2 =  @(p,x) 0+((1-0-p(1))./(1 + exp(-(x-p(2))/p(3)))); % logistic
mdl2=fitnlm(SNR,y2,sigmf2,[0.25 0 1.5]);
thresh = interp1(sigmf2(mdl2.Coefficients.Estimate,fitsnrs),fitsnrs,...
    max(sigmf2(mdl2.Coefficients.Estimate,fitsnrs))/2,'linear','extrap');


%plot the figure
figure; 
errorbar(all_stimSNRs,nanmean(pupil_dia_pct),pupil_pct_sem,'ko'); hold on
plot(fitsnrs,sigmf2(mdl2.Coefficients.Estimate,fitsnrs),'m');
line([-40,thresh],[sigmf2(mdl2.Coefficients.Estimate,thresh),sigmf2(mdl2.Coefficients.Estimate,thresh)]);
line([thresh,thresh],[0,sigmf2(mdl2.Coefficients.Estimate,thresh)]);
xlabel('SNR (dB)'); ylabel('Fraction deviant trials with sig. response')
title('Psychometric curve fit to dilation reliability')
disp(['R2_3 = ' num2str(mdl2.Rsquared.Ordinary) ' Threshold = ' num2str(thresh) ])
