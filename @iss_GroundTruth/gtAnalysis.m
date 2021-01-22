%% GadGroundTruthNewTransform 

%%
%Method = 'Pixel'; %OMP or Pixel
Method = 'OMP';
pf = o.CallMethodPrefix(Method);
        
if strcmpi('OMP',Method)
    %Get primary and secondary sets
    [~,SortedCoefs]=sort(o.ompCoefs(:,1:73)','descend');
    SortedCoefs = SortedCoefs';
    PrimarySet = o.ompSpotCodeNo==SortedCoefs(:,1);
    SecondarySet = o.ompSpotCodeNo==SortedCoefs(:,2);
    %o.ompNeighbThresh = 13;
    %o.ompIntensityThresh = 700;
elseif strcmpi('Pixel',Method)
    %Get primary and secondary sets
    PrimarySet = o.pxSpotScore>0;
    SecondarySet = o.pxSpotScore==0;
    %o.pScoreThresh = 50;
    %o.pScoreThresh2 = 10;
    %o.pIntensityThresh = 200;
    %o.pLogProbThresh2 = 0;
    %o.pIntensityThresh2 = 50;
end
QualOK = quality_threshold(o,Method);
fprintf('Total Primary Spots: \t\t\t\t%d\n',sum(QualOK&PrimarySet));
fprintf('Total Primary or Secondary Spots: \t\t%d\n',sum(QualOK&(PrimarySet|SecondarySet)));
fprintf('Total Spots: \t\t\t\t\t%d\n',sum(QualOK));

for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
        pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
        fprintf('There are %d %s peak spots\n', sum(o.gtTruePositiveSet{r,b}),...
            o.GeneNames{o.gtGeneNo(r,b)});
        fprintf('of which, we can achieve %d\n',sum(pfTruePosSet));
        fprintf('False positive set has %d spots.\n',sum(pfFalsePosSet));

        %print excel data
        fprintf('Total Primary True Positives: \t\t\t%d\n',sum(QualOK&PrimarySet&pfTruePosSet));
        fprintf('Total Primary or Secondary True Positives: \t%d\n',...
            sum(QualOK&(PrimarySet|SecondarySet)&pfTruePosSet));
        fprintf('Total True Positives: \t\t\t\t%d\n',sum(QualOK&pfTruePosSet));
        fprintf('Total Primary False Positives: \t\t\t%d\n',sum(QualOK&PrimarySet&pfFalsePosSet));
        fprintf('Total Primary or Secondary False Positives: \t%d\n',...
            sum(QualOK&(PrimarySet|SecondarySet)&pfFalsePosSet));
        fprintf('Total False Positives: \t\t\t\t%d\n',sum(QualOK&pfFalsePosSet));
        fprintf('\n');
    end
end

%% Save data
%Make Tables
MakeTables = false;
if MakeTables
    FileLocation = {''};
    Method = {''};
    Intensity_Method = {''};
    ScoreScale = 0;
    pScoreThresh = 0;
    pIntensityThresh = 0;
    pScoreThresh2 = 0;
    pIntensityThresh2 = 0;
    pLogProbThresh2 = 0;
    pQualThresh1 = 0;
    pQualParam1 = 0;
    pQualThresh2 = 0;
    pQualParam2 = 0;
    ompNeighbThresh = 0;
    ompIntensityThresh = 0;
    ompScoreThresh = 0;
    nPrimarySpots = 0;
    nPrimarySecondarySpots = 0;
    nTotalSpots = 0;
    Combined_Score = 0;
    Combined_nTP = 0;
    Combined_nTP_Missed = 0;
    Combined_nFP = 0;
    TruePosData.Summary = table(FileLocation,Method,Intensity_Method,ScoreScale,pScoreThresh,pIntensityThresh,...
        pScoreThresh2,pIntensityThresh2,pLogProbThresh2,pQualThresh1,pQualParam1,pQualThresh2,...
        pQualParam2,ompNeighbThresh,ompIntensityThresh,ompScoreThresh,nPrimarySpots,nPrimarySecondarySpots,...
        nTotalSpots,Combined_Score,Combined_nTP,Combined_nTP_Missed,Combined_nFP);
    
    TP_Max = 0;
    TP_Primary = 0;
    TP_PrimarySecondary = 0;
    TP_Total = 0;
    FP_Max = 0;
    FP_Primary = 0;
    FP_PrimarySecondary = 0;
    FP_Total = 0;
    for r=o.gtRounds
        for b=o.UseChannels
            if o.gtGeneNo(r,b)==0; continue; end
            TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b)))) = ...
                table(TP_Max,TP_Primary,TP_PrimarySecondary,TP_Total,...
                FP_Max,FP_Primary,FP_PrimarySecondary,FP_Total);
        end
    end
end
%% Add data to tables
% Run after the above so have QualOK and Primary.

i = size(TruePosData.Summary,1)+1;      %INDEX OF DATA TO BE ADDED
TruePosData.Summary(i,'FileLocation') = ...
    {fullfile(o.OutputDirectory, 'oCall_spots_OMP_Unbled')};
%Method = 'Pixel: pLogThresh, ProbMethod=2';
%Method = 'Pixel: pQualThresh, ProbMethod=2';
%Method = 'OMP: BledCodes';
Method = 'OMP: UnBledCodes';
%Method = 'Spatial';
Intensity_Method = 'Mean Unbled - Mean Not Unbled';
%Intensity_Method = 'Mean Unbled';
%Intensity_Method = 'Median Unbled';
%Intensity_Method = 'Mean Unbled - Mean Not Unbled, Z-Scored';

TruePosData.Summary(i,'Method') = {Method};
TruePosData.Summary(i,'Intensity_Method') = {Intensity_Method};
for k=1:size(TruePosData.Summary,2)
    try
        TruePosData.Summary(i,k) = {nan};
    end
end

if strcmpi(Method,'Pixel: pLogThresh') || strcmpi(Method,'Pixel: pLogThresh, ProbMethod=2')
    if o.ProbMethod==1
        TruePosData.Summary(i,'ScoreScale') = {o.ScoreScale};
    end
    TruePosData.Summary(i,'pScoreThresh') = {o.pScoreThresh};
    TruePosData.Summary(i,'pIntensityThresh') = {o.pIntensityThresh};
    TruePosData.Summary(i,'pScoreThresh2') = {o.pScoreThresh2};
    TruePosData.Summary(i,'pIntensityThresh2') = {o.pIntensityThresh2};
    TruePosData.Summary(i,'pLogProbThresh2') = {o.pLogProbThresh2};
elseif strcmpi(Method,'Pixel: pQualThresh') || strcmpi(Method,'Pixel: pQualThresh, ProbMethod=2')
    TruePosData.Summary(i,'Intensity_Method') = {'N/A'};
    if o.ProbMethod==1
        TruePosData.Summary(i,'ScoreScale') = {o.ScoreScale};
    end
    TruePosData.Summary(i,'pQualThresh1') = {o.pQualThresh1};
    TruePosData.Summary(i,'pQualParam1') = {o.pQualParam1};
    TruePosData.Summary(i,'pQualThresh2') = {o.pQualThresh2};
    TruePosData.Summary(i,'pQualParam2') = {o.pQualParam2};
elseif strcmpi(Method,'OMP: BledCodes') || strcmpi(Method,'OMP: UnBledCodes')
    TruePosData.Summary(i,'ompNeighbThresh') = {o.ompNeighbThresh};
    TruePosData.Summary(i,'ompIntensityThresh') = {o.ompIntensityThresh};
    TruePosData.Summary(i,'ompScoreThresh') = {o.ompScoreThresh};
end
TruePosData.Summary(i,'nPrimarySpots') = {sum(QualOK&PrimarySet)};
TruePosData.Summary(i,'nPrimarySecondarySpots') = {sum(QualOK&(PrimarySet|SecondarySet))};
TruePosData.Summary(i,'nTotalSpots') = {sum(QualOK)};
Combined_Score = 0;
Combined_nTP = 0;
Combined_nTP_Missed = 0;
Combined_nFP = 0;
for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'nTP') = ...
            {sum(o.gtTruePositiveSet{r,b})};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_Max') = ...
            {sum(pfTruePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_Primary') = ...
            {sum(QualOK&PrimarySet&pfTruePosSet)}; 
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_PrimarySecondary') = ...
            {sum(QualOK&(PrimarySet|SecondarySet)&pfTruePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_Total') = ...
            {sum(QualOK&pfTruePosSet)};
        
        pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_Max') = ...
            {sum(pfFalsePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_Primary') = ...
            {sum(QualOK&PrimarySet&pfFalsePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_PrimarySecondary') = ...
            {sum(QualOK&(PrimarySet|SecondarySet)&pfFalsePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_Total') = ...
            {sum(QualOK&pfFalsePosSet)};
        
        Combined_nTP = Combined_nTP+sum(QualOK&pfTruePosSet);
        Combined_nTP_Missed = Combined_nTP_Missed+sum(~QualOK&pfTruePosSet)+...
            sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet);
        Combined_nFP = Combined_nFP+sum(QualOK&pfFalsePosSet);
        %Combined_Score = Combined_Score+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
        %    5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);    
        %Min Score is best
        Combined_Score = Combined_Score+sum(~QualOK&pfTruePosSet)+...
            sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
    end
end
TruePosData.Summary(i,'Combined_Score') = {Combined_Score};
TruePosData.Summary(i,'Combined_nTP') = {Combined_nTP};
TruePosData.Summary(i,'Combined_nTP_Missed') = {Combined_nTP_Missed};
TruePosData.Summary(i,'Combined_nFP') = {Combined_nFP};

%save(fullfile(o.OutputDirectory, 'GroundTruth_Data'), 'TruePosData', '-v7.3');
%% Find Best Params - PixelBased
%QualThresh method
QualThresh1 = 3:1:10;
QualParam1 = 12:1:17;
QualThresh2 = -5:1:5;
QualParam2 = 0;
ScoreImage = zeros(length(QualThresh1),length(QualParam1),length(QualThresh2),length(QualParam2));
for i=1:length(QualThresh1)
    o.pQualThresh1 = QualThresh1(i);
    for j=1:length(QualParam1)
        o.pQualParam1 = QualParam1(j);
        for k=1:length(QualThresh2)
            o.pQualThresh2 = QualThresh2(k);
            for k2=1:length(QualParam2)
                o.pQualParam2 = QualParam2(k2);
                QualOK = quality_threshold(o,'Pixel');
                for r=o.gtRounds
                    for b=o.UseChannels
                        if o.gtGeneNo(r,b)==0; continue; end
                        pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
                        pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
                        %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
                        %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
                        ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(~QualOK&pfTruePosSet)+...
                            sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
                    end
                end
            end
        end
    end
end
%figure; imagesc(ScoreImage);
[a,b] = min(ScoreImage(:));
a
[a,b,c,d]=ind2sub(size(ScoreImage),b);
o.pQualThresh1=QualThresh1(a);
o.pQualParam1=QualParam1(b);
o.pQualThresh2=QualThresh2(c);
o.pQualParam2=QualParam2(d);
QualThresh1(a)
QualParam1(b)
QualThresh2(c)
QualParam2(d)

% %pLogThresh method
% ScoreThresh = 3;
% IntensityThresh = -550:10:-450;
% o.pScoreThresh2 = -0.000001;
% IntensityThresh2 = -inf;
% LogProbThresh2 = -5:1:5;
% % ScoreThresh = 60:1:68;
% % IntensityThresh = 0:10:180;
% % o.pScoreThresh2 = -0.000001;
% % IntensityThresh2 = -inf;
% % LogProbThresh2 = -13:1:-5;
% ScoreImage = zeros(length(ScoreThresh),length(IntensityThresh),length(IntensityThresh2),length(LogProbThresh2));
% for i=1:length(ScoreThresh)
%     o.pScoreThresh = ScoreThresh(i);
%     for j=1:length(IntensityThresh)
%         o.pIntensityThresh = IntensityThresh(j);
%         for k=1:length(IntensityThresh2)
%             o.pIntensityThresh2 = IntensityThresh2(k);
%             for k2=1:length(LogProbThresh2)
%                 o.pLogProbThresh2 = LogProbThresh2(k2);
%                 QualOK = quality_threshold(o,'Pixel');
%                 for r=o.gtRounds
%                     for b=o.UseChannels
%                         if o.gtGeneNo(r,b)==0; continue; end
%                         pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
%                         pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
%                         %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
%                         %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
%                         ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(~QualOK&pfTruePosSet)+...
%                             sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
%                     end
%                 end
%             end
%         end
%     end
% end
% %figure; imagesc(ScoreImage);
% [a,b] = min(ScoreImage(:));
% a
% [a,b,c,d]=ind2sub(size(ScoreImage),b);
% o.pScoreThresh=ScoreThresh(a);
% o.pIntensityThresh=IntensityThresh(b);
% o.pIntensityThresh2=IntensityThresh2(c);
% o.pLogProbThresh2=LogProbThresh2(d);
% ScoreThresh(a)
% IntensityThresh(b)
% IntensityThresh2(c)
% LogProbThresh2(d)

%% Testing OMP thresholds
NeighbThresh = 6:1:8;
IntensityThresh = 1050:10:1150;
ScoreThresh = 12:0.1:14;

ScoreImage = zeros(length(NeighbThresh),length(IntensityThresh),length(ScoreThresh));

for n=1:length(NeighbThresh)
    o.ompNeighbThresh=NeighbThresh(n);
    for i=1:length(IntensityThresh)
        o.ompIntensityThresh=IntensityThresh(i);
        for s=1:length(ScoreThresh)
            o.ompScoreThresh=ScoreThresh(s);
            QualOK = quality_threshold(o,'OMP');
            for r=o.gtRounds
                for b=o.UseChannels
                    if o.gtGeneNo(r,b)==0; continue; end
                    pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
                    pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
                    %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
                    %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
                    ScoreImage(n,i,s) = ScoreImage(n,i,s)+sum(~QualOK&pfTruePosSet)+...
                        sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
                end
            end
        end
    end
end
[a,b] = min(ScoreImage(:));
a
[a,b,c]=ind2sub(size(ScoreImage),b);
o.ompNeighbThresh=NeighbThresh(a);
o.ompIntensityThresh=IntensityThresh(b);
o.ompScoreThresh = ScoreThresh(c);
NeighbThresh(a)
IntensityThresh(b)
ScoreThresh(c)