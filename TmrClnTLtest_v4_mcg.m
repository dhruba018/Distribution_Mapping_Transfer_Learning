clc;    clear;      %close all

DATAPATH = sprintf('C:\\Users\\%s\\', getenv('username'));
DATAPATH = {[DATAPATH, 'Dropbox (Texas Tech)\#DataStorageSRD\METABRIC_cBioPortal\brca_metabric\'];
                            [DATAPATH, 'Dropbox (Texas Tech)\#DataStorageSRD\CCLE_cBioPortal\cellline_ccle_broad\'];
                            [DATAPATH, 'Google Drive\Data CCLE and GDSC\GDSC\Version 7\']};

FILENAME = {'data_expression_metabric.txt'; 'data_expression_median_ccle.txt';...
                        'gdsc_v7_gene_expression_processed_SRD_09_27_2018.xlsx'};
EmptyDataFlag = {'.', 'NA', 'N/A'};

%% Data Extraction...
%%% METABRIC...
dataTableGX1 = readtable([DATAPATH{1}, FILENAME{1}], 'TreatAsEmpty', EmptyDataFlag);
GX1 = dataTableGX1{:, 3:end}';                            genes1 = dataTableGX1(:, 1:2);
samples1 = dataTableGX1.Properties.VariableNames(3:end)';
[~, ind11] = sort(lower(genes1.Hugo_Symbol));     [~, ind12] = sort(lower(samples1));
genes1 = genes1(ind11, :);                                     GX1 = GX1(:, ind11);
samples1 = samples1(ind12);                                  GX1 = GX1(ind12, :);

%%% CCLE (vNew)...
dataTableGX2 = readtable([DATAPATH{2}, FILENAME{2}], 'TreatAsEmpty', EmptyDataFlag);
dataGX2 = dataTableGX2{:, 3:end}';                      genes2 = dataTableGX2(:, 1:2);
CL_list2 = dataTableGX2.Properties.VariableNames(3:end)';
CL2 = cell(numel(CL_list2), 2);                             % CL table
for k = 1:numel(CL_list2)
    C = split(CL_list2(k), '_');
    CL2(k, :) = [C(1), join(C(2:end), ' ')];
end
CL2 = array2table(CL2, 'variablenames', {'CL_Name', 'Organ'});
[~, ind2] = sort(lower(CL2.Organ));                      CL2 = CL2(ind2, :);        dataGX2 = dataGX2(ind2, :);

% BRCA data...
brCLidx2 = strcmpi(CL2.Organ, 'breast');
brCL2 = CL2{brCLidx2, 1};                                   GX2 = dataGX2(brCLidx2, :);
[~, ind21] = sort(lower(genes2.Hugo_Symbol));    [~, ind22] = sort(lower(brCL2));
genes2 = genes2(ind21, :);                                    GX2 = GX2(:, ind21);
brCL2 = brCL2(ind22);                                         GX2 = GX2(ind22, :);

% All CL data...
CL2tot = CL2.CL_Name;                                        GX2tot = dataGX2(:, ind21);

%%% GDSC (v7)...
dataTableGX3 = readtable([DATAPATH{3}, FILENAME{3}], 'TreatAsEmpty', EmptyDataFlag);
dataGX3 = dataTableGX3{:, 2:end}';                    genes3 = dataTableGX3{:, 1};
CL_list3 = dataTableGX3.Properties.VariableNames(2:end)';
CL3 = cell(numel(CL_list3), 2);                             % CL table
for k = 1:numel(CL_list3)
    C = split(CL_list3(k), '_');
    CL3(k, :) = [C(1), join(C(2:end), ' ')];
end
CL3 = array2table(CL3, 'variablenames', {'CL_Name', 'Organ'});
[~, ind3] = sort(lower(CL3.Organ));                      CL3 = CL3(ind3, :);         dataGX3 = dataGX3(ind3, :);

% BRCA data...
brCLidx3 = strcmpi(CL3.Organ, 'breast');
brCL3 = CL3{brCLidx3, 1};                                   GX3 = dataGX3(brCLidx3, :);
[~, ind31] = sort(lower(genes3));                         [~, ind32] = sort(lower(brCL3));
genes3 = genes3(ind31, :);                                    GX3 = GX3(:, ind31);
brCL3 = brCL3(ind32);                                         GX3 = GX3(ind32, :);

% All CL data...
CL3tot = CL3.CL_Name;                                        GX3tot = dataGX3(:, ind31);


%% Matched data...
matGN = [ ];        [~, matGN(:, 1), matGN(:, 2), matGN(:, 3)] =... 
                                        intersectn(lower(genes1.Hugo_Symbol), lower(genes2.Hugo_Symbol), lower(genes3), 'stable');
% samples1, brCL2, brCL3
geneSet = genes1(matGN(:, 1), :);                         GXset1 = GX1(:, matGN(:, 1));
GXset2tot = GX2tot(:, matGN(:, 2));                     GXset2 = GX2(:, matGN(:, 2));
GXset3tot = GX3tot(:, matGN(:, 3));                     GXset3 = GX3(:, matGN(:, 3));

% Missing value imputation...
Knn = 10;
if sum(isnan(GXset1(:))) > 0,            GXset1 = knnimpute(GXset1', Knn)';                  end
if sum(isnan(GXset2(:))) > 0,            GXset2 = knnimpute(GXset2', Knn)';                 end
if sum(isnan(GXset3(:))) > 0,            GXset3 = knnimpute(GXset3', Knn)';                 end
if sum(isnan(GXset2tot(:))) > 0,        GXset2tot = knnimpute(GXset2tot', Knn)';        end
if sum(isnan(GXset3tot(:))) > 0,        GXset3tot = knnimpute(GXset3tot', Knn)';        end

% % Standardization (z-score)...
% GXset1 = zscore(GXset1, [ ], 1);
% GXset2 = zscore(GXset2, [ ], 1);                           GXset2tot = zscore(GXset2tot, [ ], 1);
% GXset3 = zscore(GXset3, [ ], 1);                           GXset3tot = zscore(GXset3tot, [ ], 1);

% Finding biomarkers in data...
bmChoice = 2;
switch bmChoice
    case 1
        biomarkerList = importdata('Biomarker_list.csv');
        [~, ind] = intersect(lower(geneSet.Hugo_Symbol), lower(biomarkerList), 'stable');
    case 2
        biomarkerList = readtable('Biomarker_list2.csv');
        [~, ind] = intersect(lower(geneSet.Hugo_Symbol), lower(biomarkerList.Symbol), 'stable');
end
bmIdx = false(numel(geneSet.Hugo_Symbol), 1);       bmIdx(ind) = true;
biomarkers = geneSet.Hugo_Symbol(bmIdx);             XgeneSet = geneSet(~bmIdx, :);

% Biomarker expression as response...
Xdata1 = GXset1(:, ~bmIdx);                                Ydata1 = GXset1(:, bmIdx);
Xdata2 = GXset2(:, ~bmIdx);                               Ydata2 = GXset2(:, bmIdx);
Xdata3 = GXset3(:, ~bmIdx);                               Ydata3 = GXset3(:, bmIdx);
Xdata2tot = GXset2tot(:, ~bmIdx);                      Ydata2tot = GXset2tot(:, bmIdx);
Xdata3tot = GXset3tot(:, ~bmIdx);                      Ydata3tot = GXset3tot(:, bmIdx);
p = sum(~bmIdx);                                                 q = sum(bmIdx);
n1 = numel(samples1);
n2 = numel(brCL2);                                              n2tot = numel(CL2tot);
n3 = numel(brCL3);                                              n3tot = numel(CL3tot);

% Standardization (z-score)...
GXset11 = zscore(GXset1, [ ], 1);
GXset22 = zscore(GXset2, [ ], 1);                         GXset22tot = zscore(GXset2tot, [ ], 1);
GXset33 = zscore(GXset3, [ ], 1);                         GXset33tot = zscore(GXset3tot, [ ], 1);
Xdata11 = zscore(Xdata1, [ ], 1);                          Ydata11 = zscore(Ydata1, [ ], 1);
Xdata22 = zscore(Xdata2, [ ], 1);                         Ydata22 = zscore(Ydata2, [ ], 1);
Xdata33 = zscore(Xdata3, [ ], 1);                         Ydata33 = zscore(Ydata3, [ ], 1);
Xdata22tot = zscore(Xdata2tot, [ ], 1);                Ydata22tot = zscore(Ydata2tot, [ ], 1);
Xdata33tot = zscore(Xdata3tot, [ ], 1);                Ydata33tot = zscore(Ydata3tot, [ ], 1);


%%% Distributions...
numBin = 50;

% Biomarkers...
hb = [figure(2), figure(3)];            nPts = 1e3;            clf(figure(hb(1)));       clf(figure(hb(2)));
mbox = [1:q/2, 3*q/2+1:2*q; q/2+1:q, 2*q+1:5*q/2; q+1:3*q/2, 5*q/2+1:3*q];
for k = 1 : q
    [Ypdf1, pts1] = ksdensity(Ydata11(:, k), 'NumPoints', nPts);
    [Ypdf2, pts2] = ksdensity(Ydata22(:, k), 'NumPoints', nPts);
    [Ypdf3, pts3] = ksdensity(Ydata33(:, k), 'NumPoints', nPts);
    [Ypdf2tot, pts2tot] = ksdensity(Ydata22tot(:, k), 'NumPoints', nPts);
    [Ypdf3tot, pts3tot] = ksdensity(Ydata33tot(:, k), 'NumPoints', nPts);
    
    figure(hb(1))               % BRCA
    subplot(6,q/2,mbox(1,k)),     histogram(Ydata11(:, k), numBin, 'normalization', 'pdf')
    hold on,                            plot(pts1, Ypdf1, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0.5,0,0}METABRIC'])
    subplot(6,q/2,mbox(2,k)),     histogram(Ydata22(:, k), numBin, 'normalization', 'pdf')
    hold on,                             plot(pts2, Ypdf2, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0.5,0}CCLE'])
    subplot(6,q/2,mbox(3,k)),     histogram(Ydata33(:, k), numBin, 'normalization', 'pdf')
    hold on,                             plot(pts3, Ypdf3, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0,0.5}GDSC'])
   
    figure(hb(2))               % All CLs
    subplot(6,q/2,mbox(1,k)),     histogram(Ydata11(:, k), numBin, 'normalization', 'pdf')
    hold on,                            plot(pts1, Ypdf1, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0.5,0,0}METABRIC'])
    subplot(6,q/2,mbox(2,k)),     histogram(Ydata22tot(:, k), numBin, 'normalization', 'pdf')
    hold on,                             plot(pts2tot, Ypdf2tot, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0.5,0}CCLE'])
    subplot(6,q/2,mbox(3,k)),     histogram(Ydata33tot(:, k), numBin, 'normalization', 'pdf')
    hold on,                             plot(pts3tot, Ypdf3tot, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0,0.5}GDSC'])
end
defn = ['Row 1, 3 - {\color[rgb]{0,0,0.8}METABRIC}, Row 2, 4 - {\color[rgb]{0,0,0.8}CCLE}, ',... 
                    'Row 3, 6 - {\color[rgb]{0,0,0.8}GDSC}'];
figure(hb(1)),           suptitle([{'\bfBiomarker Distributions for {\color[rgb]{0.5,0,0}BRCA} CLs'}; {defn}])
figure(hb(2)),          suptitle([{'\bfBiomarker Distributions for {\color[rgb]{0.5,0,0}All} CLs'}; {defn}])


%% Housekeeping Gene Data...
HKGtab = readtable('HKG_list.csv');
nHKG = HKGtab.nHKG;             nHKG(cellfun(@isempty, nHKG)) = [ ];
tHKG = HKGtab.tHKG;              tHKG(cellfun(@isempty, tHKG)) = [ ];
HKGlist = {nHKG, tHKG, [nHKG; tHKG]};
[~, nHKGidx] = intersect(lower(XgeneSet.Hugo_Symbol), lower(nHKG), 'stable');
[~, tHKGidx] = intersect(lower(XgeneSet.Hugo_Symbol), lower(tHKG), 'stable');
HKGidx = {nHKGidx, tHKGidx};        pHKG = [numel(HKGidx{1}), numel(HKGidx{2})];    % #HKGs
HKGset = {XgeneSet.Hugo_Symbol(HKGidx{1}), XgeneSet.Hugo_Symbol(HKGidx{2})};

% Distributions...
numBin = 50;        nPts = 1e3;
defn = {'nHKG', 'tHKG', ['Row 1 - {\color[rgb]{0,0,0.8}METABRIC}, ',...
                'Row 2 - {\color[rgb]{0,0,0.8}CCLE: BRCA},  Row 3 - {\color[rgb]{0,0,0.8}GDSC: BRCA}']};
hh = [figure(4), figure(5)];            clf(figure(hh(1)));       clf(figure(hh(2)));
for i = 1 : 2
    figure(hh(i))
    for k = 1 : pHKG(i)
        [Xpdf1, pts1]  = ksdensity(GXset11(:, HKGidx{i}(k)), 'NumPoints', nPts);
        [Xpdf2, pts2] = ksdensity(GXset22(:, HKGidx{i}(k)), 'NumPoints', nPts);
        [Xpdf3, pts3] = ksdensity(GXset33(:, HKGidx{i}(k)), 'NumPoints', nPts);
        
        subplot(3,7,k),                  histogram(GXset11(:, HKGidx{i}(k)), numBin, 'normalization', 'pdf')
        hold on,                            plot(pts1, Xpdf1, 'r-.', 'linewidth', 2),         hold off
        title([HKGset{i}{k}, ' - \color[rgb]{0.5,0,0}METABRIC'])
        subplot(3,7,k+7),               histogram(GXset22(:, HKGidx{i}(k)), numBin, 'normalization', 'pdf')
        hold on,                             plot(pts2, Xpdf2, 'r-.', 'linewidth', 2),       hold off
        title([HKGset{i}{k}, ' - \color[rgb]{0,0.5,0}CCLE: BRCA'])
        subplot(3,7,k+14),              histogram(GXset33(:, HKGidx{i}(k)), numBin, 'normalization', 'pdf')
        hold on,                             plot(pts3, Xpdf3, 'r-.', 'linewidth', 2),       hold off
        title([HKGset{i}{k}, ' - \color[rgb]{0,0,0.5}GDSC: BRCA'])
    end
    suptitle([{['\bfHousekeeping Gene Distributions - {\color[rgb]{0.5,0,0}', defn{i}, '}']; defn{3}}])
end


%% Pairwise Dependency...
seed = 0;       rng(seed)

% Non-biomarker data...
repSwitch = 2;
switch repSwitch
    case 1
        fprintf(1, 'Chosen non-biomarkers = '),        fprintf(2, 'Random regular genes\n')
        gnIdx = randperm(p, q);                               % Random representatives, # = q
    case 2
        hkSwitch = 3;       hkChoice = {'nHKG', 'tHKG', 'All'};
        fprintf(1, 'Chosen non-biomarkers = '),        fprintf(2, '%s\n', hkChoice{hkSwitch})
        housekeepingList = HKGlist{hkSwitch};
        [~, gnIdx] = intersect(lower(XgeneSet.Hugo_Symbol), lower(housekeepingList), 'stable');
end
repGenes = XgeneSet.Hugo_Symbol(gnIdx);        pRep = numel(repGenes);
XdataRep1 = Xdata1(:, gnIdx);
XdataRep2 = Xdata2(:, gnIdx);                            XdataRep2tot = Xdata2tot(:, gnIdx);
XdataRep3 = Xdata3(:, gnIdx);                            XdataRep3tot = Xdata3tot(:, gnIdx);

% Gene-pair combinations...
bmPairs = nchoosek(1:q, 2);       nPairsBM = size(bmPairs, 1);
bmPairList = strrep(cellstr([num2str(bmPairs(:, 1)), repmat('-', [nPairsBM, 1]), num2str(bmPairs(:, 2))]), ' ', '');
gnPairs = nchoosek(1:pRep, 2);  nPairsGN = size(gnPairs, 1);
gnPairList = strrep(cellstr([num2str(gnPairs(:, 1)), repmat('-', [nPairsGN, 1]), num2str(gnPairs(:, 2))]), ' ', '');


%%% Pairwise Dependency - Correlation...
% Biomarkers...
ycorrmat1 = corr(Ydata1, 'type', 'spearman');                  ycorr1 = 1 - squareform(1 - ycorrmat1, 'tovector')';
ycorrmat2 = corr(Ydata2, 'type', 'spearman');                 ycorr2 = 1 - squareform(1 - ycorrmat2, 'tovector')';
ycorrmat3 = corr(Ydata3, 'type', 'spearman');                 ycorr3 = 1 - squareform(1 - ycorrmat3, 'tovector')';
ycorrmat2tot = corr(Ydata2tot, 'type', 'spearman');         ycorr2tot = 1 - squareform(1 - ycorrmat2tot, 'tovector')';
ycorrmat3tot = corr(Ydata3tot, 'type', 'spearman');         ycorr3tot = 1 - squareform(1 - ycorrmat3tot, 'tovector')';
ycorr_struct = [corr(ycorr1, ycorr2, 'type', 'pearson'), corr(ycorr1, ycorr2tot, 'type', 'pearson');...
                            corr(ycorr1, ycorr3, 'type', 'pearson'), corr(ycorr1, ycorr3tot, 'type', 'pearson');...
                            corr(ycorr2, ycorr3, 'type', 'pearson'), corr(ycorr2tot, ycorr3tot, 'type', 'pearson')];

% Non-biomarker genes...
xcorrmat1 = corr(XdataRep1, 'type', 'spearman');               xcorr1 = 1 - squareform(1 - xcorrmat1, 'tovector')';
xcorrmat2 = corr(XdataRep2, 'type', 'spearman');              xcorr2 = 1 - squareform(1 - xcorrmat2, 'tovector')';
xcorrmat3 = corr(XdataRep3, 'type', 'spearman');              xcorr3 = 1 - squareform(1 - xcorrmat3, 'tovector')';
xcorrmat2tot = corr(XdataRep2tot, 'type', 'spearman');      xcorr2tot = 1 - squareform(1 - xcorrmat2tot, 'tovector')';
xcorrmat3tot = corr(XdataRep3tot, 'type', 'spearman');      xcorr3tot = 1 - squareform(1 - xcorrmat3tot, 'tovector')';
xcorr_struct = [corr(xcorr1, xcorr2, 'type', 'pearson'), corr(xcorr1, xcorr2tot, 'type', 'pearson');...
                            corr(xcorr1, xcorr3, 'type', 'pearson'), corr(xcorr1, xcorr3tot, 'type', 'pearson');...
                            corr(xcorr2, xcorr3, 'type', 'pearson'), corr(xcorr2tot, xcorr3tot, 'type', 'pearson')];

% Results...
fprintf(1, 'Correlation betn Gene-pair corr. vectors for Biomarkers & ')
switch repSwitch
    case 1
        fprintf(1, 'Regular Genes:\n')
    case 2
        fprintf(1, 'Housekeeping Genes (%s):\n', hkChoice{hkSwitch})
end
rNames = {'r_MC - BM'; 'r_MG - BM'; 'r_CG - BM'; 'r_MC'; 'r_MG'; 'r_CG'};
corr_struct =  array2table(round([ycorr_struct; xcorr_struct], 4), 'rownames', rNames,...
                                            'variablenames', {'BRCA', 'ALL'});
disp(corr_struct)

% Correlation plots...
switch repSwitch
    case 1
        hc = figure(7);                                             % Regular genes
    case 2
        switch hkSwitch
            case 1
                hc = figure(8);                                     % nHKG
            case 2
                hc = figure(9);                                     % tHKG
            case 3
                hc = figure(10);                                   % All HKGs
        end
end
defn = {'Regular Genes'; 'Housekeeping Genes'; ' (nHKG)'; ' (tHKG)'; ''};
defn2 = {'Row 1 - For {\color[rgb]{0,0,0.8}BRCA} CLs, Row 2 - For {\color[rgb]{0,0,0.8}All} CLs'};

% Biomarkers...
figure(6),                        clf
subplot(241),                   plot(1:nPairsBM, ycorr1, 'rs--', 1:nPairsBM, ycorr2, 'gs--', 1:nPairsBM, ycorr3, 'bs--')
hold on,                           line([0, nPairsBM+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--')
xlim([0, nPairsBM+1]),      xticksNoted = floor(linspace(1, nPairsBM, 12));
xticks(xticksNoted),        xticklabels(bmPairList(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),   legend({'\bfMETABRIC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}BRCA}'}; {['({\it\color[rgb]{0.5,0,0}r_{MC} = ',... 
                num2str(round(ycorr_struct(1,1), 2)), ', r_{MG} = ', num2str(round(ycorr_struct(2,1), 2)),... 
                ', r_{CG} = ', num2str(round(ycorr_struct(3,1), 2)), '})']}])
subplot(242),                  scatter(ycorr1, ycorr2, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),      ylabel('\bfCCLE')
title('{\color[rgb]{0,0,0.8}METABRIC vs. CCLE} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(243),                  scatter(ycorr1, ycorr3, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),      ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}METABRIC vs. GDSC} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(244),                  scatter(ycorr2, ycorr3, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),                ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}CCLE vs. GDSC} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(245),                  plot(1:nPairsBM, ycorr1, 'rs--', 1:nPairsBM, ycorr2tot, 'gs--', 1:nPairsBM, ycorr3tot, 'bs--')
hold on,                          line([0, nPairsBM+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--')
xlim([0, nPairsBM+1]),     xticksNoted = floor(linspace(1, nPairsBM, 12));
xticks(xticksNoted),       xticklabels(bmPairList(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),   legend({'\bfMETABRIC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}ALL}'}; {['({\it\color[rgb]{0.5,0,0}r_{MC} = ',... 
                num2str(round(ycorr_struct(1,2), 2)), ', r_{MG} = ', num2str(round(ycorr_struct(2,2), 2)),... 
                ', r_{CG} = ', num2str(round(ycorr_struct(3,2), 2)), '})']}])
subplot(246),                  scatter(ycorr1, ycorr2tot, 24, 'bo', 'filled'),                 box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),      ylabel('\bfCCLE')
title('{\color[rgb]{0,0,0.8}METABRIC vs. CCLE} [{\color[rgb]{0.5,0.4,0}ALL}]')
subplot(247),                  scatter(ycorr1, ycorr3tot, 24, 'bo', 'filled'),                 box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),      ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}METABRIC vs. GDSC} [{\color[rgb]{0.5,0.4,0}ALL}]')
subplot(248),                  scatter(ycorr2tot, ycorr3tot, 24, 'bo', 'filled'),            box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),                ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}CCLE vs. GDSC} [{\color[rgb]{0.5,0.4,0}ALL}]')
suptitle([{'\bfComparison of Gene-pair Correlations for {\color[rgb]{0.5,0,0}Biomarkers}'}; defn2])

% Non-biomarker genes...
figure(hc),                      clf
subplot(241),                   plot(1:nPairsGN, xcorr1, 'rs--', 1:nPairsGN, xcorr2, 'gs--', 1:nPairsGN, xcorr3, 'bs--')
hold on,                           line([0, nPairsGN+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--'),      ylim([-0.4, 0.6])
xlim([0, nPairsGN+1]),      xticksNoted = floor(linspace(1, nPairsGN, 12));
xticks(xticksNoted),        xticklabels(gnPairList(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),   legend({'\bfMETABRIC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}BRCA}'}; {['({\it\color[rgb]{0.5,0,0}r_{MC} = ',... 
                num2str(round(xcorr_struct(1,1), 2)), ', r_{MG} = ', num2str(round(xcorr_struct(2,1), 2)),... 
                ', r_{CG} = ', num2str(round(xcorr_struct(3,1), 2)), '})']}])
subplot(242),                  scatter(xcorr1, xcorr2, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),       ylabel('\bfCCLE')
title('{\color[rgb]{0,0,0.8}METABRIC vs. CCLE} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(243),                  scatter(xcorr1, xcorr3, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),       ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}METABRIC vs. GDSC} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(244),                  scatter(xcorr2, xcorr3, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),       ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}CCLE vs. GDSC} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(245),                  plot(1:nPairsGN, xcorr1, 'rs--', 1:nPairsGN, xcorr2tot, 'gs--', 1:nPairsGN, xcorr3tot, 'bs--')
hold on,                          line([0, nPairsGN+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--'),      ylim([-0.4, 0.6])
xlim([0, nPairsGN+1]),     xticksNoted = floor(linspace(1, nPairsGN, 12));
xticks(xticksNoted),       xticklabels(gnPairList(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),   legend({'\bfMETABRIC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}ALL}'}; {['({\it\color[rgb]{0.5,0,0}r_{MC} = ',... 
                num2str(round(xcorr_struct(1,2), 2)), ', r_{MG} = ', num2str(round(xcorr_struct(2,2), 2)),... 
                ', r_{CG} = ', num2str(round(xcorr_struct(3,2), 2)), '})']}])
subplot(246),                  scatter(xcorr1, xcorr2tot, 24, 'bo', 'filled'),                 box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),       ylabel('\bfCCLE')
title('{\color[rgb]{0,0,0.8}METABRIC vs. CCLE} [{\color[rgb]{0.5,0.4,0}ALL}]')
subplot(247),                  scatter(xcorr1, xcorr3tot, 24, 'bo', 'filled'),                 box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),       ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}METABRIC vs. GDSC} [{\color[rgb]{0.5,0.4,0}ALL}]')
subplot(248),                  scatter(xcorr2tot, xcorr3tot, 24, 'bo', 'filled'),            box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),                ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}CCLE vs. GDSC} [{\color[rgb]{0.5,0.4,0}ALL}]')
if repSwitch == 1
    suptitle([{['\bfComparison of Gene-pair Correlations for {\color[rgb]{0.5,0,0}',  defn{repSwitch}, '}']}; defn2])
elseif repSwitch == 2
    suptitle([{['\bfComparison of Gene-pair Correlations for {\color[rgb]{0.5,0,0}', defn{repSwitch},... 
                        defn{hkSwitch+2}, '}']}; defn2])
end


%% Correlation Heatmaps...
% Unaltered Corr. Matrix Heatmaps...
ColMap = jet;
figure(11),           clf
subplot(221),       heatmap(ycorrmat1, 'XData', biomarkers, 'YData', biomarkers, 'Colormap', ColMap);
                           colorbar off,        title('METABRIC')
subplot(222),       heatmap(ycorrmat2, 'XData', biomarkers, 'YData', biomarkers, 'Colormap', ColMap);
                           title('CCLE: BRCA')
subplot(224),       heatmap(ycorrmat3, 'XData', biomarkers, 'YData', biomarkers, 'Colormap', ColMap);
                           title('GDSC: BRCA')
suptitle('\bfBiomarker Correlation Structure')

% Hierarchical clustering...
ydistmat1 = squareform(pdist(Ydata1', 'euclidean'));                 yTree1 = linkage(ydistmat1, 'average');
ydistmat2 = squareform(pdist(Ydata2', 'euclidean'));                yTree2 = linkage(ydistmat2, 'average');
ydistmat3 = squareform(pdist(Ydata3', 'euclidean'));                yTree3 = linkage(ydistmat3, 'average');
ydistmat2tot = squareform(pdist(Ydata2tot', 'euclidean'));        yTree2tot = linkage(ydistmat2tot, 'average');
ydistmat3tot = squareform(pdist(Ydata3tot', 'euclidean'));        yTree3tot = linkage(ydistmat3tot, 'average');
nCls = 3;                                                              % #clusters
IdxCls1 = cluster(yTree1, 'maxclust', nCls);
IdxCls2 = cluster(yTree2, 'maxclust', nCls);        IdxCls2tot = cluster(yTree2tot, 'maxclust', nCls);
IdxCls3 = cluster(yTree3, 'maxclust', nCls);        IdxCls3tot = cluster(yTree3tot, 'maxclust', nCls);
w = 1.000;              TreeCut = zeros(3, 2);
TreeCut(1, 1) = w * yTree1((q - nCls + 1), 3);
TreeCut(2, :) = [w * yTree2((q - nCls + 1), 3), w * yTree2tot((q - nCls + 1), 3)];
TreeCut(3, :) = [w * yTree3((q - nCls + 1), 3), w * yTree3tot((q - nCls + 1), 3)];
% TreeCut

% Cluster tree plot w/ cut-off... 
figure(12),          clf
subplot(231),      ht1 = dendrogram(yTree1);         set(ht1, 'color', 'k'),        box on
hold on,              line([0, q+1], [TreeCut(1, 1), TreeCut(1, 1)], 'linestyle', '--', 'linewidth', 1.2, 'color', 'r')
xticklabels(biomarkers(xticks)),        xtickangle(45),         ylabel('\bfLinkage Distance'),     title(['METABRIC'])
subplot(232),      ht2 = dendrogram(yTree2);       set(ht2, 'color', 'k'),        box on
hold on,              line([0, q+1], [TreeCut(2, 1), TreeCut(2, 1)], 'linestyle', '--', 'linewidth', 1.2, 'color', 'r')
xticklabels(biomarkers(xticks)),        xtickangle(45),         ylabel('\bfLinkage Distance'),     title(['CCLE: BRCA'])
subplot(233),      ht3 = dendrogram(yTree3);       set(ht3, 'color', 'k'),        box on
hold on,              line([0, q+1], [TreeCut(3, 1), TreeCut(3, 1)], 'linestyle', '--', 'linewidth', 1.2, 'color', 'r')
xticklabels(biomarkers(xticks)),        xtickangle(45),         ylabel('\bfLinkage Distance'),     title(['GDSC: BRCA'])

% Sorted clusters w/ corresponding corr. matrices...
bmSorted = cell(q, 5);                                                       [~, ind1] = sort(IdxCls1);
[~, ind2] = sort(IdxCls2);                                                  [~, ind2tot] = sort(IdxCls2tot);
[~, ind3] = sort(IdxCls3);                                                  [~, ind3tot] = sort(IdxCls3tot);
bmSorted(:, 1) = biomarkers(ind2);                                    YdataSorted1 = Ydata1(:, ind2);
bmSorted(:, 2) = biomarkers(ind2);                                    YdataSorted2 = Ydata2(:, ind2);
bmSorted(:, 3) = biomarkers(ind2);                                    YdataSorted3 = Ydata3(:, ind2);
bmSorted(:, 4) = biomarkers(ind2);                                    YdataSorted2tot = Ydata2tot(:, ind2);
bmSorted(:, 5) = biomarkers(ind2);                                    YdataSorted3tot = Ydata3tot(:, ind2);
ys_cmat1 = corr(YdataSorted1, 'type', 'spearman');           ys_corr1 = 1 - squareform(1 - ys_cmat1, 'tovector')';
ys_cmat2 = corr(YdataSorted2, 'type', 'spearman');          ys_corr2 = 1 - squareform(1 - ys_cmat2, 'tovector')';
ys_cmat3 = corr(YdataSorted3, 'type', 'spearman');          ys_corr3 = 1 - squareform(1 - ys_cmat3, 'tovector')';
ys_cmat2tot = corr(YdataSorted2tot, 'type', 'spearman');  ys_corr2tot = 1 - squareform(1 - ys_cmat2tot, 'tovector')';
ys_cmat3tot = corr(YdataSorted3tot, 'type', 'spearman');  ys_corr3tot = 1 - squareform(1 - ys_cmat3tot, 'tovector')';
ys_corr_struct = [corr(ys_corr1, ys_corr2, 'type', 'pearson'), corr(ys_corr1, ys_corr2tot, 'type', 'pearson');...
                            corr(ys_corr1, ys_corr3, 'type', 'pearson'), corr(ys_corr1, ys_corr3tot, 'type', 'pearson');...
                            corr(ys_corr2, ys_corr3, 'type', 'pearson'), corr(ys_corr2tot, ys_corr3tot, 'type', 'pearson')];

% New heatmaps...
subplot(234),       heatmap(ys_cmat1, 'XData', bmSorted(:, 1), 'YData', bmSorted(:, 1), 'Colormap', ColMap);
                            colorbar off,        title('METABRIC')
subplot(235),       heatmap(ys_cmat2, 'XData', bmSorted(:, 2), 'YData', bmSorted(:, 2), 'Colormap', ColMap);
                            colorbar off,        title('CCLE: BRCA')
subplot(236),       heatmap(ys_cmat3, 'XData', bmSorted(:, 3), 'YData', bmSorted(:, 3), 'Colormap', ColMap);
                            title('GDSC: BRCA')
defn = [sprintf('Row 1 - {\\color[rgb]{0,0,0.8}Hierarchical Clustering (K = %d)}, ', nCls),...
                            'Row 2 - {\color[rgb]{0,0,0.8}Correlation Heatmap}'];
suptitle([{'\bfBiomarker Correlation Structure'}; {defn}])


%% Distribution Smoothing w/ Spline Fitting & Kernel Smoothing Density Estimate...
%   Ydata1, Ydata2, Ydata2tot

% Cubic spline fitting...
samp0 = {linspace(0, 1, n1)', linspace(0, 1, n2)', linspace(0, 1, n3)', linspace(0, 1, n2tot)', linspace(0, 1, n3tot)'};
nPts = n1;                                                             samp1 = linspace(0, 1, nPts)';
YdataFit11 = zeros(nPts, q);
YdataFit22 = zeros(nPts, q);                                YdataFit22tot = zeros(nPts, q);
YdataFit33 = zeros(nPts, q);                                YdataFit33tot = zeros(nPts, q);
for k = 1 : q
    YdataFit11(:, k) = spline(samp0{1}, sort(Ydata11(:, k)), samp1);
    YdataFit22(:, k) = spline(samp0{2}, sort(Ydata22(:, k)), samp1);
    YdataFit33(:, k) = spline(samp0{3}, sort(Ydata33(:, k)), samp1);
    YdataFit22tot(:, k) = spline(samp0{4}, sort(Ydata22tot(:, k)), samp1);
    YdataFit33tot(:, k) = spline(samp0{5}, sort(Ydata33tot(:, k)), samp1);
end

% Smoothed distributions w/ kernel density overlayed...
hb = [figure(14), figure(15)];            clf(figure(hb(1)));       clf(figure(hb(2)));
mbox = [1:q/2, 3*q/2+1:2*q; q/2+1:q, 2*q+1:5*q/2; q+1:3*q/2, 5*q/2+1:3*q];          nPts = 1e3;
for k = 1 : q
    [Ypdf1, pts1] = ksdensity(Ydata11(:, k), 'NumPoints', nPts);
    [Ypdf2, pts2] = ksdensity(Ydata22(:, k), 'NumPoints', nPts);
    [Ypdf3, pts3] = ksdensity(Ydata33(:, k), 'NumPoints', nPts);
    [Ypdf2tot, pts2tot] = ksdensity(Ydata22tot(:, k), 'NumPoints', nPts);
    [Ypdf3tot, pts3tot] = ksdensity(Ydata33tot(:, k), 'NumPoints', nPts);
    
    figure(hb(1))                                                   % BRCA
    subplot(6,q/2,mbox(1,k)),     histogram(YdataFit11(:, k), numBin, 'normalization', 'pdf')
    hold on,                                plot(pts1, Ypdf1, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0.5,0,0}METABRIC'])
    subplot(6,q/2,mbox(2,k)),     histogram(YdataFit22(:, k), numBin, 'normalization', 'pdf')
    hold on,                                plot(pts2, Ypdf2, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0.5,0}CCLE'])
    subplot(6,q/2,mbox(3,k)),     histogram(YdataFit33(:, k), numBin, 'normalization', 'pdf')
    hold on,                                plot(pts3, Ypdf3, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0,0.5}GDSC'])
   
    figure(hb(2))                                                   % All CLs
    subplot(6,q/2,mbox(1,k)),     histogram(YdataFit11(:, k), numBin, 'normalization', 'pdf')
    hold on,                                plot(pts1, Ypdf1, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0.5,0,0}METABRIC'])
    subplot(6,q/2,mbox(2,k)),     histogram(YdataFit22tot(:, k), numBin, 'normalization', 'pdf')
    hold on,                                plot(pts2tot, Ypdf2tot, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0.5,0}CCLE'])
    subplot(6,q/2,mbox(3,k)),     histogram(YdataFit33tot(:, k), numBin, 'normalization', 'pdf')
    hold on,                                plot(pts3tot, Ypdf3tot, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0,0.5}GDSC'])
end
defn = ['Row 1, 3 - {\color[rgb]{0,0,0.8}METABRIC}, Row 2, 4 - {\color[rgb]{0,0,0.8}CCLE}, ',... 
                    'Row 3, 6 - {\color[rgb]{0,0,0.8}GDSC}'];
figure(hb(1)),      suptitle([{['\bf{\color[rgb]{0.5,0,0}Spline-fitted} Biomarker Distributions ',... 
                                                    'for {\color[rgb]{0.5, 0.4, 0}BRCA} CLs']}; {defn}])
figure(hb(2)),      suptitle([{['\bf{\color[rgb]{0.5,0,0}Spline-fitted} Biomarker Distributions ',... 
                                                    'for {\color[rgb]{0.5, 0.4, 0}All} CLs']}; {defn}])


%% Feature selection - TIME-CONSUMING LOOP!!!...
% RELIEFF gene ranking...
K = 10;                           rank1 = zeros(p, q);
rank2 = zeros(p, q);        rank2tot = zeros(p, q);
rank3 = zeros(p, q);        rank3tot = zeros(p, q);
tic
for k = 1 : q
    rank1(:, k) = relieff(Xdata1, Ydata1(:, k), K, 'method', 'regression');
    rank2(:, k) = relieff(Xdata2, Ydata2(:, k), K, 'method', 'regression');
    rank3(:, k) = relieff(Xdata3, Ydata3(:, k), K, 'method', 'regression');
    rank2tot(:, k) = relieff(Xdata2tot, Ydata2tot(:, k), K, 'method', 'regression');
    rank3tot(:, k) = relieff(Xdata3tot, Ydata3tot(:, k), K, 'method', 'regression');
end
FEATFILENAME = sprintf('relieff_ranked_genes_MCG_v2_bm%d.mat', bmChoice);
save(FEATFILENAME, 'rank1', 'rank2', 'rank3', 'rank2tot', 'rank3tot')
toc


%% Pairwise Dependency - Top Genes...
% Load RELIEFF ranking...
FEATFILENAME = sprintf('relieff_ranked_genes_MCG_v2_bm%d.mat', bmChoice);
load(FEATFILENAME, 'rank1', 'rank2', 'rank3', 'rank2tot', 'rank3tot')

% Finding common genes [Gene names => XgeneSet.Hugo_Symbol]...
nFeat = 10e3;           rankImp = cell(q, 1);       rankImpTot = cell(q, 1);        geneImp = cell(q, 2);
fprintf(1, 'No. of top features used = '),        fprintf(2, '%d\n', nFeat)
fprintf(1, 'Intersection of common features betn all %d biomarkers:\n', q)
cRank = (1 : p)';       cRankTot = (1 : p)';
for k = 1 : q
    rankImp{k} = intersectn(rank1(1:nFeat, k), rank2(1:nFeat, k), rank3(1:nFeat, k), 'sorted');
    rankImpTot{k} = intersectn(rank1(1:nFeat, k), rank2tot(1:nFeat, k), rank3tot(1:nFeat, k), 'sorted');
    geneImp(k, :) = {XgeneSet.Hugo_Symbol(rankImp{k}), XgeneSet.Hugo_Symbol(rankImpTot{k})};
    cRank = intersectn(cRank, rankImp{k}, 'sorted');
    cRankTot = intersectn(cRankTot, rankImpTot{k}, 'sorted');
    % fprintf(1, 'Size: BRCA = %d x %d,  ALL = %d x %d\n', size(cRank), size(cRankTot))
end
sizeImp = array2table([cellfun(@length, rankImp), cellfun(@length, rankImpTot)],...
                                    'variablenames', {'BRCA', 'ALL'}, 'rownames', biomarkers);

cGeneImp = XgeneSet.Hugo_Symbol(cRank);                  pFS = numel(cGeneImp);
cGeneImpTot = XgeneSet.Hugo_Symbol(cRankTot);        pFStot = numel(cGeneImpTot);
cRank2 = intersect(cRank, cRankTot, 'sorted');     % Overall common genes for both BRCA & All CLs
cGeneCL = XgeneSet.Hugo_Symbol(cRank2);                   pFS2 = numel(cGeneCL);

% Print results...
% fprintf(1, 'No. of common features betn tumor & cell lines...\n'),      disp(sizeImp)
fprintf(1, 'No. of common important genes:\n')
disp(array2table([pFS, pFStot, pFS2]', 'rownames', {'BRCA'; 'ALL'; 'Intersection'}, 'variablenames', {'nCL'}))

% Common gene data [BRCA, ALL]...
XdataFS1 = Xdata1(:, cRank);           XdataFS1tot = Xdata1(:, cRankTot);
XdataFS2 = Xdata2(:, cRank);          XdataFS2tot = Xdata2tot(:, cRankTot);
XdataFS3 = Xdata3(:, cRank);          XdataFS3tot = Xdata3tot(:, cRankTot);

% % Pairwise correlation...
% Pairs...
gnPairsFS = nchoosek(1 : pFS, 2);               nPairsFS = size(gnPairsFS, 1);
gnPairListFS = strrep(cellstr([num2str(gnPairsFS(:, 1)), repmat('-', [nPairsFS, 1]), num2str(gnPairsFS(:, 2))]), ' ', '');
gnPairsFStot = nchoosek(1 : pFStot, 2);       nPairsFStot = size(gnPairsFStot, 1);
gnPairListFStot = strrep(cellstr([num2str(gnPairsFStot(:, 1)), repmat('-', [nPairsFStot, 1]),...
                                                    num2str(gnPairsFStot(:, 2))]), ' ', '');

% Correlation...
xcorrmatFS1 = corr(XdataFS1, 'type', 'spearman');         xcorrFS1 = 1 - squareform(1 - xcorrmatFS1, 'tovector')';
xcorrmatFS2 = corr(XdataFS2, 'type', 'spearman');        xcorrFS2 = 1 - squareform(1 - xcorrmatFS2, 'tovector')';
xcorrmatFS3 = corr(XdataFS3, 'type', 'spearman');        xcorrFS3 = 1 - squareform(1 - xcorrmatFS3, 'tovector')';
xcorrmatFS1tot = corr(XdataFS1tot, 'type', 'spearman');
xcorrFS1tot = 1 - squareform(1 - xcorrmatFS1tot, 'tovector')';
xcorrmatFS2tot = corr(XdataFS2tot, 'type', 'spearman');
xcorrFS2tot = 1 - squareform(1 - xcorrmatFS2tot, 'tovector')';
xcorrmatFS3tot = corr(XdataFS3tot, 'type', 'spearman');
xcorrFS3tot = 1 - squareform(1 - xcorrmatFS3tot, 'tovector')';

xcorrFS_struct = [corr(xcorrFS1, xcorrFS2, 'type', 'pearson'), corr(xcorrFS1tot, xcorrFS2tot, 'type', 'pearson');...
                            corr(xcorrFS1, xcorrFS3, 'type', 'pearson'), corr(xcorrFS1tot, xcorrFS3tot, 'type', 'pearson');...
                            corr(xcorrFS2, xcorrFS3, 'type', 'pearson'), corr(xcorrFS2tot, xcorrFS3tot, 'type', 'pearson')];

% Display results...
fprintf(1, 'Correlation betn Gene-pair corr. vectors:\n')
disp(array2table([pFS, pFStot; round(xcorrFS_struct, 4)], 'variablenames', {'BRCA', 'ALL'},....
                                                        'rownames', {'nGenesTop'; 'r_MC'; 'r_MG'; 'r_CG'}))

figure(16)
subplot(241),                   plot(1:nPairsFS, xcorrFS1, 'rs--', 1:nPairsFS, xcorrFS2, 'gs--',... 
                                                1:nPairsFS, xcorrFS3, 'bs--')
hold on,                           line([0, nPairsFS+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--'),     % ylim([-0.4, 0.6])
xlim([0, nPairsFS+1]),      xticksNoted = floor(linspace(1, nPairsFS, 12));
xticks(xticksNoted),        xticklabels(gnPairListFS(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),   legend({'\bfMETABRIC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}BRCA}'}; {['({\it\color[rgb]{0.5,0,0}r_{MC} = ',... 
                num2str(round(xcorrFS_struct(1,1), 2)), ', r_{MG} = ', num2str(round(xcorrFS_struct(2,1), 2)),... 
                ', r_{CG} = ', num2str(round(xcorrFS_struct(3,1), 2)), '})']}])
subplot(242),                  scatter(xcorrFS1, xcorrFS2, 24, 'bo', 'filled'),              box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),        ylabel('\bfCCLE')
title('{\color[rgb]{0,0,0.8}METABRIC vs. CCLE} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(243),                  scatter(xcorrFS1, xcorrFS3, 24, 'bo', 'filled'),              box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),        ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}METABRIC vs. GDSC} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(244),                  scatter(xcorrFS2, xcorrFS3, 24, 'bo', 'filled'),             box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),                 ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}CCLE vs. GDSC} [{\color[rgb]{0.5,0.4,0}BRCA}]')
subplot(245),                  plot(1:nPairsFStot, xcorrFS1tot, 'rs--', 1:nPairsFStot, xcorrFS2tot, 'gs--',... 
                                                1:nPairsFStot, xcorrFS3tot, 'bs--')
hold on,                          line([0, nPairsFStot+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--'),     % ylim([-0.4, 0.6])
xlim([0, nPairsFStot+1]), xticksNoted = floor(linspace(1, nPairsFStot, 12));
xticks(xticksNoted),       xticklabels(gnPairListFStot(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),   legend({'\bfMETABRIC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}ALL}'}; {['({\it\color[rgb]{0.5,0,0}r_{MC} = ',... 
                num2str(round(xcorrFS_struct(1,2), 2)), ', r_{MG} = ', num2str(round(xcorrFS_struct(2,2), 2)),... 
                ', r_{CG} = ', num2str(round(xcorrFS_struct(3,2), 2)), '})']}])
subplot(246),                  scatter(xcorrFS1tot, xcorrFS2tot, 24, 'bo', 'filled'),     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),        ylabel('\bfCCLE')
title('{\color[rgb]{0,0,0.8}METABRIC vs. CCLE} [{\color[rgb]{0.5,0.4,0}ALL}]')
subplot(247),                  scatter(xcorrFS1tot, xcorrFS3tot, 24, 'bo', 'filled'),     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfMETABRIC'),        ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}METABRIC vs. GDSC} [{\color[rgb]{0.5,0.4,0}ALL}]')
subplot(248),                  scatter(xcorrFS2tot, xcorrFS3tot, 24, 'bo', 'filled'),     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),                 ylabel('\bfGDSC')
title('{\color[rgb]{0,0,0.8}CCLE vs. GDSC} [{\color[rgb]{0.5,0.4,0}ALL}]')
defn = {'Row 1 - For {\color[rgb]{0,0,0.8}BRCA} CLs, Row 2 - For {\color[rgb]{0,0,0.8}All} CLs'};
suptitle([{'\bfComparison of Gene-pair Correlations for {\color[rgb]{0.5,0,0}Top Genes}'}; defn])


%% Additional analyses...
norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));

nPts = 1e3;                                   pts = linspace(0, 1, nPts)';
bmPairs = nchoosek(1:q, 2);          nPairsBM = size(bmPairs, 1);         bmPairsNameList = join(biomarkers(bmPairs), ' - ');
Zstates = {'s00', 's01', 's10', 's11'};                    % 4-state discretization
Zpmf1 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Zpmf2 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Zpmf3 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Th = [0.25, 0.75];                                                % Distribution tail quantiles
asc_met = 'NMI';                                                % Association metric
YpairLo_asc = NaN(nPairsBM, 3);         YpairUp_asc = NaN(nPairsBM, 3);
for b = 1 : nPairsBM
    Ypair1 = norm01(Ydata1(:, bmPairs(b, :)));        Ypair_pdf1 = [ksdensity(Ypair1(:, 1), pts), ksdensity(Ypair1(:, 2), pts)];
    Ypair2 = norm01(Ydata2(:, bmPairs(b, :)));       Ypair_pdf2 = [ksdensity(Ypair2(:, 1), pts), ksdensity(Ypair2(:, 2), pts)];
    Ypair3 = norm01(Ydata3(:, bmPairs(b, :)));       Ypair_pdf3 = [ksdensity(Ypair3(:, 1), pts), ksdensity(Ypair3(:, 2), pts)];
    
    % Otsu's thresholding...
    ot1 = [otsuthresh(Ypair_pdf1(:, 1)), otsuthresh(Ypair_pdf1(:, 2))];
    ot2 = [otsuthresh(Ypair_pdf2(:, 1)), otsuthresh(Ypair_pdf2(:, 2))];
    ot3 = [otsuthresh(Ypair_pdf3(:, 1)), otsuthresh(Ypair_pdf3(:, 2))];
    
    % Discretization...
    z1 = ~(Ypair1 <= ot1);        Zpmf1{b, :} = [mean(sum(z1 == [0, 0], 2) == 2), mean(sum(z1 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z1 == [1, 0], 2) == 2), mean(sum(z1 == [1, 1], 2) == 2)];
    z2 = ~(Ypair2 <= ot2);      Zpmf2{b, :} = [mean(sum(z2 == [0, 0], 2) == 2), mean(sum(z2 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z2 == [1, 0], 2) == 2), mean(sum(z2 == [1, 1], 2) == 2)];
    z3 = ~(Ypair3 <= ot3);      Zpmf3{b, :} = [mean(sum(z3 == [0, 0], 2) == 2), mean(sum(z3 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z3 == [1, 0], 2) == 2), mean(sum(z3 == [1, 1], 2) == 2)];
    
    % Distribution tail threshold quantiles...
    qt1 = quantile(Ypair1, Th);         qt2 = quantile(Ypair2, Th);         qt3 = quantile(Ypair3, Th);
    
    % Lower & upper dist. tails...
    YpairLo1 = Ypair1((sum(Ypair1 < qt1(1, :), 2) == 2), :);          YpairHi1 = Ypair1((sum(Ypair1 > qt1(2, :), 2) == 2), :);
    YpairLo2 = Ypair2((sum(Ypair2 < qt2(1, :), 2) == 2), :);        YpairHi2 = Ypair2((sum(Ypair2 > qt2(2, :), 2) == 2), :);
    YpairLo3 = Ypair3((sum(Ypair3 < qt3(1, :), 2) == 2), :);        YpairHi3 = Ypair3((sum(Ypair3 > qt3(2, :), 2) == 2), :);
    
    % Association betn. pairwise distribution tails...
    switch asc_met
        case 'corr'
            % Correlation betn. tails...
            if numel(YpairLo1) > 2,     YpairLo_asc(b, 1) = corr(YpairLo1(:, 1), YpairLo1(:, 2), 'type', 'spearman');       end
            if numel(YpairLo2) > 2,     YpairLo_asc(b, 2) = corr(YpairLo2(:, 1), YpairLo2(:, 2), 'type', 'spearman');     end
            if numel(YpairLo3) > 2,     YpairLo_asc(b, 3) = corr(YpairLo3(:, 1), YpairLo3(:, 2), 'type', 'spearman');     end
            if numel(YpairHi1) > 2,      YpairUp_asc(b, 1) = corr(YpairHi1(:, 1), YpairHi1(:, 2), 'type', 'spearman');     end
            if numel(YpairHi2) > 2,     YpairUp_asc(b, 2) = corr(YpairHi2(:, 1), YpairHi2(:, 2), 'type', 'spearman');    end
            if numel(YpairHi3) > 2,     YpairUp_asc(b, 3) = corr(YpairHi3(:, 1), YpairHi3(:, 2), 'type', 'spearman');    end
            defn = 'Correlation';
        case 'NMI'
            % Normalized Mutual Info betn. tails...
            %%% [MI, Idx] = MutualInformation(X, opt, scale, dtype, C)
            if numel(YpairLo1) > 2,     YpairLo_asc(b, 1) = MutualInformation(YpairLo1, 'single', 'corr');       end
            if numel(YpairLo2) > 2,     YpairLo_asc(b, 2) = MutualInformation(YpairLo2, 'single', 'corr');      end
            if numel(YpairLo3) > 2,     YpairLo_asc(b, 3) = MutualInformation(YpairLo3, 'single', 'corr');      end
            if numel(YpairHi1) > 2,      YpairUp_asc(b, 1) = MutualInformation(YpairHi1, 'single', 'corr');      end
            if numel(YpairHi2) > 2,     YpairUp_asc(b, 2) = MutualInformation(YpairHi2, 'single', 'corr');      end
            if numel(YpairHi3) > 2,     YpairUp_asc(b, 3) = MutualInformation(YpairHi3, 'single', 'corr');      end
            defn = 'Normalized MI';
    end
    
    % Printing results...
    if ~mod(b, 10)
        pmfDiff = array2table(round([Zpmf1{b, :} - Zpmf2{b, :}; Zpmf1{b, :} - Zpmf3{b, :}; Zpmf2{b, :} - Zpmf3{b, :}], 4),...
                                                    'RowNames', {'M-C'; 'M-G'; 'C-G'}, 'VariableNames', Zstates);
        fprintf(1, '\nBiomarker pair [%d] = ', b),        fprintf(2, '%s\n', bmPairsNameList{b})
        fprintf(1, 'Error betn discrete PMFs = \n'),	disp(pmfDiff)
        ascTail = array2table(round([YpairLo_asc(b, :); YpairUp_asc(b, :)], 4), 'RowNames', {'Lower tail'; 'Upper tail'},...
                                                        'VariableNames', {'METABRIC', 'CCLE', 'GDSC'});
        fprintf(1, '%s between lower & upper distribution tails = \n', asc_met),        disp(ascTail)
    end
end

% %%
% Visualization...
numBin = 50;
bmPairList = strrep(cellstr([num2str(bmPairs(:, 1)), repmat('-', [nPairsBM, 1]), num2str(bmPairs(:, 2))]), ' ', '');

figure(17),         clf
ax1 = axes('Position', [0.05, 0.14, 0.15, 0.6]);       ax2 = axes('Position', [0.21, 0.75, 0.65, 0.2]);
ax3 = axes('Position', [0.21, 0.14, 0.65, 0.6], 'YAxisLocation', 'right');
histogram(ax1, Ypair1(:, 2), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax1, [-90, 90]),	hold on
histogram(ax2, Ypair1(:, 1), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax2, [0, 90]),      hold on
scatter(ax3, Ypair1(:, 1), Ypair1(:, 2), 18, 'filled', 'o'),       box on
hold on,       line(ax3, [0, 1], [ot1(2), ot1(2)], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt1(1, 2), qt1(1, 2)], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt1(1, 1), qt1(1, 1)], [0, 1], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [ot1(1), ot1(1)], [0, 1], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt1(2, 2), qt1(2, 2)], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt1(2, 1), qt1(2, 1)], [0, 1], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
legend('\bf\itData', '\bfTh_{OTSU}', '\bfTh_{LO}', '\bfTh_{UP}')
xlabel(sprintf('\\bf%s', biomarkers{bmPairs(b, 1)})),      ylabel(sprintf('\\bf%s', biomarkers{bmPairs(b, 2)}))
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}METABRIC}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')
% saveas(figure(17), ['.\Figures__MB_CCLE_GDSC_Analysis\v2\Biomarker_pair_scatters\',... 
%             sprintf('Biomarker_pair_scatter_%s_%s_MB.png', biomarkers{bmPairs(b, 1)}, biomarkers{bmPairs(b, 2)})])

figure(18),         clf
ax1 = axes('Position', [0.05, 0.14, 0.15, 0.6]);       ax2 = axes('Position', [0.21, 0.75, 0.65, 0.2]);
ax3 = axes('Position', [0.21, 0.14, 0.65, 0.6], 'YAxisLocation', 'right');
histogram(ax1, Ypair2(:, 2), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax1, [-90, 90]),	hold on
histogram(ax2, Ypair2(:, 1), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax2, [0, 90]),      hold on
scatter(ax3, Ypair2(:, 1), Ypair2(:, 2), 18, 'filled', 'o'),       box on
hold on,       line(ax3, [0, 1], [ot2(2), ot2(2)], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt2(1, 2), qt2(1, 2)], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt2(1, 1), qt2(1, 1)], [0, 1], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [ot2(1), ot2(1)], [0, 1], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt2(2, 2), qt2(2, 2)], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt2(2, 1), qt2(2, 1)], [0, 1], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
legend('\bf\itData', '\bfTh_{OTSU}', '\bfTh_{LO}', '\bfTh_{UP}')
xlabel(sprintf('\\bf%s', biomarkers{bmPairs(b, 1)})),      ylabel(sprintf('\\bf%s', biomarkers{bmPairs(b, 2)}))
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}CCLE}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')
% saveas(figure(18), ['.\Figures__MB_CCLE_GDSC_Analysis\v2\Biomarker_pair_scatters\',... 
%             sprintf('Biomarker_pair_scatter_%s_%s_CCLE.png', biomarkers{bmPairs(b, 1)}, biomarkers{bmPairs(b, 2)})])

figure(19),         clf
ax1 = axes('Position', [0.05, 0.14, 0.15, 0.6]);       ax2 = axes('Position', [0.21, 0.75, 0.65, 0.2]);
ax3 = axes('Position', [0.21, 0.14, 0.65, 0.6], 'YAxisLocation', 'right');
histogram(ax1, Ypair3(:, 2), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax1, [-90, 90]),	hold on
histogram(ax2, Ypair3(:, 1), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax2, [0, 90]),      hold on
scatter(ax3, Ypair3(:, 1), Ypair3(:, 2), 18, 'filled', 'o'),       box on
hold on,       line(ax3, [0, 1], [ot3(2), ot3(2)], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt3(1, 2), qt3(1, 2)], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt3(1, 1), qt3(1, 1)], [0, 1], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [ot3(1), ot3(1)], [0, 1], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt3(2, 2), qt3(2, 2)], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt3(2, 1), qt3(2, 1)], [0, 1], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
legend('\bf\itData', '\bfTh_{OTSU}', '\bfTh_{LO}', '\bfTh_{UP}')
xlabel(sprintf('\\bf%s', biomarkers{bmPairs(b, 1)})),      ylabel(sprintf('\\bf%s', biomarkers{bmPairs(b, 2)}))
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}GDSC}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')
% saveas(figure(19), ['.\Figures__MB_CCLE_GDSC_Analysis\v2\Biomarker_pair_scatters\',... 
%             sprintf('Biomarker_pair_scatter_%s_%s_GDSC.png', biomarkers{bmPairs(b, 1)}, biomarkers{bmPairs(b, 2)})])
% end

% %
% Linear model fit...
model12 = {fitlm(Zpmf1{:, 1}, Zpmf2{:, 1}), fitlm(Zpmf1{:, 2}, Zpmf2{:, 2}),...
                        fitlm(Zpmf1{:, 3}, Zpmf2{:, 3}), fitlm(Zpmf1{:, 4}, Zpmf2{:, 4})};
model13 = {fitlm(Zpmf1{:, 1}, Zpmf3{:, 1}), fitlm(Zpmf1{:, 2}, Zpmf3{:, 2}),... 
                        fitlm(Zpmf1{:, 3}, Zpmf3{:, 3}), fitlm(Zpmf1{:, 4}, Zpmf3{:, 4})};
w12 = [model12{1}.Coefficients.Estimate, model12{2}.Coefficients.Estimate,... 
                model12{3}.Coefficients.Estimate, model12{4}.Coefficients.Estimate];
w13 = [model13{1}.Coefficients.Estimate, model13{2}.Coefficients.Estimate,... 
                model13{3}.Coefficients.Estimate, model13{4}.Coefficients.Estimate];
Wm = [mean([w12(1, :); w13(1, :)], 1); mean([w12(2, :); w13(2, :)], 1)];                        % Mean coefficient
xval = linspace(0, round(max(max(Zpmf1{:, :})), 1), 100)';       Xval = [ones(size(xval)), xval];
yfit = [Xval * Wm(:, 1), Xval * Wm(:, 2), Xval * Wm(:, 3), Xval * Wm(:, 4)];                % Fitted line

% Discretized biomarker-pair PMFs...
figure(20),          clf
subplot(221),       scatter(Zpmf1{:, 1}, Zpmf2{:, 1}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 1}, Zpmf3{:, 1}, 24, 'ro', 'filled'),                    box on
hold on,               plot(xval, yfit(:, 1), 'g-'),          xlabel('\bfTumor Culture'),       ylabel('\bfCell line')
legend('\bfMETABRIC - CCLE', '\bfMETABRIC - GDSC', 'Location', 'nw'),          title('State = \color[rgb]{0,0,0.8}00')
subplot(222),       scatter(Zpmf1{:, 2}, Zpmf2{:, 2}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 2}, Zpmf3{:, 2}, 24, 'ro', 'filled'),                   box on
hold on,               plot(xval, yfit(:, 2), 'g-'),          xlabel('\bfTumor Culture'),      ylabel('\bfCell line')
legend('\bfMETABRIC - CCLE', '\bfMETABRIC - GDSC', 'Location', 'nw'),         title('State = \color[rgb]{0,0,0.8}01')
subplot(223),       scatter(Zpmf1{:, 3}, Zpmf2{:, 3}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 3}, Zpmf3{:, 3}, 24, 'ro', 'filled'),                  box on
hold on,               plot(xval, yfit(:, 3), 'g-'),          xlabel('\bfTumor Culture'),     ylabel('\bfCell line')
legend('\bfMETABRIC - CCLE', '\bfMETABRIC - GDSC', 'Location', 'nw'),        title('State = \color[rgb]{0,0,0.8}10')
subplot(224),       scatter(Zpmf1{:, 4}, Zpmf2{:, 4}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 4}, Zpmf3{:, 4}, 24, 'ro', 'filled'),                 box on
hold on,               plot(xval, yfit(:, 4), 'g-'),          xlabel('\bfTumor Culture'),    ylabel('\bfCell line')
legend('\bfMETABRIC - CCLE', '\bfMETABRIC - GDSC', 'Location', 'nw'),        title('State = \color[rgb]{0,0,0.8}11')
suptitle(['\bfComparison of {\color[rgb]{0,0,0.8}Discretized} Gene-pair Fractions for Biomarkers ',...
                    '[{\color[rgb]{0.6,0.2,0.2}BRCA}]'])
% saveas(figure(20), '.\Figures__MB_CCLE_GDSC_Analysis\v2\Biomarker_pair_discrete_fractions_BRCA.png')

% Corr. betn. tails...
figure(21),      clf
subplot(231),   scatter(YpairLo_asc(:, 1), YpairLo_asc(:, 2), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 1))+0.1, -0.1, max(YpairLo_asc(:, 2))+0.1])
xlabel('\bfMETABRIC'),         ylabel('\bfCCLE'),       title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(232),   scatter(YpairLo_asc(:, 1), YpairLo_asc(:, 3), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 1))+0.1, -0.1, max(YpairLo_asc(:, 3))+0.1])
xlabel('\bfMETABRIC'),         ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(233),   scatter(YpairLo_asc(:, 2), YpairLo_asc(:, 3), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 2))+0.1, -0.1, max(YpairLo_asc(:, 3))+0.1])
xlabel('\bfCCLE'),                  ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(234),   scatter(YpairUp_asc(:, 1), YpairUp_asc(:, 2), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 1))+0.1, -0.1, max(YpairUp_asc(:, 2))+0.1])
xlabel('\bfMETABRIC'),         ylabel('\bfCCLE'),       title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
subplot(235),   scatter(YpairUp_asc(:, 1), YpairUp_asc(:, 3), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 1))+0.1, -0.1, max(YpairUp_asc(:, 3))+0.1])
xlabel('\bfMETABRIC'),         ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
subplot(236),   scatter(YpairUp_asc(:, 2), YpairUp_asc(:, 3), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 2))+0.1, -0.1, max(YpairUp_asc(:, 3))+0.1])
xlabel('\bfCCLE'),                  ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
suptitle('\bfComparison of Gene-pair Distribution {\itTails} for Biomarkers [{\color[rgb]{0.6,0.2,0.2}BRCA}]')
% saveas(figure(21), '.\Figures__MB_CCLE_GDSC_Analysis\v2\Biomarker_pair_tail_associations_BRCA.png')

% Correlation betn. association vectors...
kLo12 = ~(isnan(YpairLo_asc(:, 1)) | isnan(YpairLo_asc(:, 2)));
kLo13 = ~(isnan(YpairLo_asc(:, 1)) | isnan(YpairLo_asc(:, 3)));
kLo23 = ~(isnan(YpairLo_asc(:, 2)) | isnan(YpairLo_asc(:, 3)));
kUp12 = ~(isnan(YpairUp_asc(:, 1)) | isnan(YpairUp_asc(:, 2)));
kUp13 = ~(isnan(YpairUp_asc(:, 1)) | isnan(YpairUp_asc(:, 3)));
kUp23 = ~(isnan(YpairUp_asc(:, 2)) | isnan(YpairUp_asc(:, 3)));
YpairTail_ascStruct = array2table(round([corr(YpairLo_asc(kLo12, 1), YpairLo_asc(kLo12, 2), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo13, 1), YpairLo_asc(kLo13, 3), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo23, 2), YpairLo_asc(kLo23, 3), 'type', 'spearman');...
                                            corr(YpairUp_asc(kUp12, 1), YpairUp_asc(kUp12, 2), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp13, 1), YpairUp_asc(kUp13, 3), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp23, 2), YpairUp_asc(kUp23, 3), 'type', 'spearman')], 4),...
                                        'rownames', {'Lower tail'; 'Upper tail'}, 'variablename', {'METABRIC', 'CCLE', 'GDSC'});
disp(YpairTail_ascStruct)


%% TRANSFER LEARNING v4...
% Xdata1, Xdata2, Ydata1, Ydata2

% % Directory info for image saving...
% SAVEPATH = sprintf('C:\\Users\\%s\\', getenv('username'));
% SAVEPATH = [SAVEPATH, 'Dropbox (Texas Tech)\Tumor_CL_TL_2018\Figures__MB_CCLE_GDSC_Analysis\v2\'];

% Load RELEIFF features...
FEATFILENAME = sprintf('relieff_ranked_genes_MCG_v2_bm%d.mat', bmChoice);
load(FEATFILENAME, 'rank1', 'rank2', 'rank3', 'rank2tot', 'rank3tot')

% Metric function definitions...
norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));       % Normalize in [0, 1]
NMAE = @(y, y_hat) mean(abs(y - y_hat)) / mean(abs(y - mean(y)));
NRMSE = @(y, y_hat) sqrt(mean((y - y_hat).^2)) / std(y, 1);

fprintf(2, 'Tumor cultures to cell line data Transfer Learning...\n')
PredApp = {'DMTL', 'DMTL_SS', 'CATL', 'BL'};
RESULTS = array2table(zeros(q+1, numel(PredApp)), 'RowNames', [biomarkers; {'Mean'}], 'VariableNames', PredApp);
RESULTS = struct('NRMSE', RESULTS, 'NMAE', RESULTS, 'SCC', RESULTS, 'numFeat', zeros(q+1, 1));
for chosenBMidx = 1 : q
    %%% Feature selection...
    fprintf(1, 'Chosen biomarker = '),       fprintf(2, '%s\n', biomarkers{chosenBMidx})
    ranks = [rank1(:, chosenBMidx), rank2(:, chosenBMidx), rank3(:, chosenBMidx)];
    
%     fprintf(1, 'Extracting common CL features...\n')                    % Extracting a common set of size m_opt
    nGN = 300;                  m_opt = 150;        gnRank = intersectn(ranks(1:nGN, 2), ranks(1:nGN, 3), 'sorted');
    m = numel(gnRank);       nI = 0;                 m0 = m;
    while m < m_opt
        nI = nI + 1;               nGN = nGN + 100;
        gnRank = intersectn(ranks(1:nGN, 2), ranks(1:nGN, 3), 'sorted');        m = numel(gnRank);
    end
%     fprintf(1, '\t#genes used = %d, #Iterations = %d\n\tInitial size = %d, Final size = %d\n', nGN, nI, m0, m)
    
    % Defining primary & secondary sets...
%     fprintf(2, 'Tumor cultures to cell line data Transfer Learning...\n')
    dsChoice = 1;                                                   % Switch for dataset choices
    switch dsChoice
        case 1                                                          % Tumor => primary, CL => secondary
            dsFlag = {'METABRIC', 'CCLE + GDSC'};
            gnList = XgeneSet.Hugo_Symbol(gnRank);         m = numel(gnList);
            X1 = Xdata1(:, gnRank);                                     X2 = [Xdata2(:, gnRank); Xdata3(:, gnRank)];
            Y1 = (Ydata1(:, chosenBMidx));                            
            Y2 = [(Ydata2(:, chosenBMidx)); (Ydata3(:, chosenBMidx))];
        case 2                                                         % CL => primary, Tumor => secondary
            dsFlag = {'CCLE + GDSC', 'METABRIC'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 1));                              m = numel(gnList);
            X1 = [Xdata2(:, ranks(1:m_opt, 1)); Xdata3(:, ranks(1:m_opt, 1))];          X2 = Xdata1(:, ranks(1:m_opt, 1));
            Y1 = [Ydata2(:, chosenBMidx); Ydata3(:, chosenBMidx)];                      Y2 = Ydata1(:, chosenBMidx);
        case 3                                                         % Tumor => primary, CCLE => secondary
            dsFlag = {'METABRIC', 'CCLE'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 2));       m = numel(gnList);
            X1 = Xdata1(:, ranks(1:m_opt, 2));                                   X2 = Xdata2(:, ranks(1:m_opt, 2));
            Y1 = Ydata1(:, chosenBMidx);                                         Y2 = Ydata2(:, chosenBMidx);
        case 4                                                         % CCLE => primary, Tumor => secondary
            dsFlag = {'CCLE', 'METABRIC'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 1));        m = numel(gnList);
            X1 = Xdata2(:, ranks(1:m_opt, 1));                                   X2 = Xdata1(:, ranks(1:m_opt, 1));
            Y1 = Ydata2(:, chosenBMidx);                                         Y2 = Ydata1(:, chosenBMidx);
        case 5                                                         % Tumor => primary, GDSC => secondary
            dsFlag = {'METABRIC', 'GDSC'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 3));       m = numel(gnList);
            X1 = Xdata1(:, ranks(1:m_opt, 3));                                   X2 = Xdata3(:, ranks(1:m_opt, 3));
            Y1 = Ydata1(:, chosenBMidx);                                         Y2 = Ydata3(:, chosenBMidx);
        case 6                                                         % GDSC => primary, Tumor => secondary
            dsFlag = {'GDSC', 'METABRIC'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 1));       m = numel(gnList);
            X1 = Xdata3(:, ranks(1:m_opt, 1));                                  X2 = Xdata1(:, ranks(1:m_opt, 1));
            Y1 = Ydata3(:, chosenBMidx);                                        Y2 = Ydata1(:, chosenBMidx);
    end
    assert(size(X1, 2) == size(X2, 2))                     % Check if #features are same
    X1 = (X1);        Y1 = norm01(Y1);         N1 = numel(Y1);
    X2 = (X2);       Y2 = norm01(Y2);        N2 = numel(Y2);
    fprintf(1, 'Primary = '),       fprintf(2, '%s,\t', dsFlag{1})
    fprintf(1, 'Secondary = '),  fprintf(2, '%s\n', dsFlag{2})
    
    % Save BM datasets for analyses...
    datapath = sprintf('C:\\Users\\%s\\Dropbox (Texas Tech)\\#DataStorageSRD\\TumorCLfeatureSelectedData',... 
                                        getenv('username'));
    if ~mod(dsChoice, 2),	 dsChar = {'_Rev', strrep(dsFlag{1}, ' + ', '_'), 'Target'};
    else,                              dsChar = {'', strrep(dsFlag{2}, ' + ', '_'), 'Source'};        end
    datadir = sprintf('BRCA_%s_%s', dsChar{3}, dsChar{2});
    filename = sprintf('Data_BM_%d_f%d_%s%s.mat', chosenBMidx, m_opt, dsChar{2}, dsChar{1});
    save(fullfile(datapath, datadir, filename), 'X1', 'X2', 'Y1', 'Y2')
    
    % Distributions...
    figure(31),                        clf
    subplot(221),                    histogram(X1, numBin, 'normalization', 'probability')
    title(sprintf('X_{primary},  dim = \\color[rgb]{0, 0, 0.8}%d x %d', size(X1)))
    subplot(222),                    histogram(Y1, numBin, 'normalization', 'probability')
    title(sprintf('y_{primary},  dim = \\color[rgb]{0, 0, 0.8}%d x %d', size(Y1)))
    subplot(223),                    histogram(X2, numBin, 'normalization', 'probability')
    title(sprintf('X_{secondary},  dim = \\color[rgb]{0, 0, 0.8}%d x %d', size(X2)))
    subplot(224),                    histogram(Y2, numBin, 'normalization', 'probability')
    title(sprintf('y_{secondary},  dim = \\color[rgb]{0, 0, 0.8}%d x %d', size(Y2)))
    defn = sprintf('Row 1 - {\\color[rgb]{0,0,0.8}%s}, Row 2 - {\\color[rgb]{0,0,0.8}%s}', dsFlag{1}, dsFlag{2});
    suptitle([{sprintf('\\bfPrimary & Secondary Data Histograms, Biomarker = {\\color[rgb]{0, 0.5, 0.4}%s}',...
                        biomarkers{chosenBMidx})}; {defn}])
    % saveas(figure(31), [SAVEPATH, sprintf('TLP_data_distribution_%s.png', biomarkers{chosenBMidx})])
    % saveas(figure(31), [SAVEPATH, sprintf('TLP_data_distribution_MC_%s.png', biomarkers{chosenBMidx})])
    
    %%% Transfer Learning predictions...
    [Y1c, X2a, CorrMatCA] = CorrAlignTransLearn(X1, X2, Y2);
    [Y1m, Y2m, Map, densityPts, CorrMatDM] = DistMatchTransLearn(X1, Y1, X2, Y2);
    
    
    %%% GX Corr. Heatmaps...
    % Distribution matching...
    CorrStructDM = [corr(CorrMatDM.Xps1D, CorrMatDM.Xp1D, 'type', 'pearson'),...
                                        corr(CorrMatDM.Xps1D, CorrMatDM.Xs1D, 'type', 'pearson')];
    figure(32),        clf;         ColMap = jet;
    subplot(221),     heatmap(CorrMatDM.Xp, 'XLabel', {'Top Genes'}, 'YLabel', {'Top Genes'}, 'Colormap', ColMap);
    colorbar off;     title(['Primary: ', dsFlag{1}])
    subplot(223),     heatmap(CorrMatDM.Xs, 'XLabel', {'Top Genes'}, 'YLabel', {'Top Genes'}, 'Colormap', ColMap);
    colorbar off;     title(['Secondary: ', dsFlag{2}])
    subplot(224),     heatmap(CorrMatDM.Xps, 'XLabel', {'Top Genes'}, 'YLabel', {'Top Genes'}, 'Colormap', ColMap);
    title('Primary Data mapped to SS')
    subplot(222),    str = {'{\bfCorrelation between...}';...
        ['    Mapped primary & primary Corr vectors = {\bf', num2str(round(CorrStructDM(1), 4)), '}'];...
        ['    Mapped primary & secondary Corr vectors = {\bf', num2str(round(CorrStructDM(2), 4)), '}']};
    dim = [0.6, 0.5, 0.5, 0.3];       ann = annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    ann.FitBoxToText = 'on';      ann.EdgeColor = 'b';      ann.LineWidth = 2;        ann.FontSize = 12;    axis off
    suptitle(['\bfCovariate Correlation Structure, Biomarker = ', biomarkers{chosenBMidx}])
    % % saveas(figure(32), ['./Figures__MB_CCLE_GDSC_Analysis/',...
    % %                                     'Covariate_heatmaps_before_after_mapping_MCG_', biomarkers{chosenBMidx}, '.png'])
    % % saveas(figure(32), ['./Figures__MB_CCLE_GDSC_Analysis/',...
    % %                                     'Covariate_heatmaps_before_after_mapping_MCG_', biomarkers{chosenBMidx}, '.fig'])
    % saveas(figure(32),...
    %     [SAVEPATH, sprintf('TLP_covariate_heatmaps_before_after_mapping_%s.png', biomarkers{chosenBMidx})])
    % saveas(figure(32),...
    %     [SAVEPATH, sprintf('TLP_covariate_heatmaps_before_after_mapping_MC_%s.png', biomarkers{chosenBMidx})])
    
    % Correlation alignment...
    CorrStructCA = [corr(CorrMatCA.Xp1D, CorrMatCA.Xsp1D, 'type', 'pearson'),...
                                        corr(CorrMatCA.Xp1D, CorrMatCA.Xs1D, 'type', 'pearson')];
    figure(33),        clf;         ColMap = jet;
    subplot(221),     heatmap(CorrMatCA.Xp, 'XLabel', {'Top Genes'}, 'YLabel', {'Top Genes'}, 'Colormap', ColMap);
    colorbar off;     title(['Primary: ', dsFlag{1}])
    subplot(223),     heatmap(CorrMatCA.Xs, 'XLabel', {'Top Genes'}, 'YLabel', {'Top Genes'}, 'Colormap', ColMap);
    colorbar off;     title(['Secondary: ', dsFlag{2}])
    subplot(224),     heatmap(CorrMatCA.Xsp, 'XLabel', {'Top Genes'}, 'YLabel', {'Top Genes'}, 'Colormap', ColMap);
    title('Secondary Data after CORAL')
    subplot(222),    str = {'{\bfCorrelation between...}';...
        ['    Aligned secondary & primary Corr vectors = {\bf', num2str(round(CorrStructCA(1), 4)), '}'];...
        ['    Aligned secondary & secondary Corr vectors = {\bf', num2str(round(CorrStructCA(2), 4)), '}']};
    dim = [0.6, 0.5, 0.5, 0.3];       ann = annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    ann.FitBoxToText = 'on';      ann.EdgeColor = 'b';      ann.LineWidth = 2;        ann.FontSize = 12;    axis off
    suptitle(['\bfCovariate Correlation Structure, Biomarker = ', biomarkers{chosenBMidx}])
    % % saveas(figure(33), ['./Figures__MB_CCLE_GDSC_Analysis/',...
    % %                                     'Covariate_heatmaps_before_after_aligning_MCG_', biomarkers{chosenBMidx}, '.png'])
    % % saveas(figure(33), ['./Figures__MB_CCLE_GDSC_Analysis/',...
    % %                                     'Covariate_heatmaps_before_after_aligning_MCG_', biomarkers{chosenBMidx}, '.fig'])
    % saveas(figure(33),...
    %     [SAVEPATH, sprintf('TLP_covariate_heatmaps_before_after_aligning_%s.png', biomarkers{chosenBMidx})])
    % saveas(figure(33),...
    %     [SAVEPATH, sprintf('TLP_covariate_heatmaps_before_after_aligning_MC_%s.png', biomarkers{chosenBMidx})])
    %%%
    
    %%% Baseline model prediction...
    TrainX = zscore(X2);        TrainY = Y2;       TestX = zscore(X1);
    rng(0);                              nTree = 200;       RF = TreeBagger(nTree, TrainX, TrainY, 'method', 'regression');
    PredY = predict(RF, TestX);                         Y1b = PredY;        %Y1b(Y1b < 0) = 0;           Y1b(Y1b > 1) = 1;
    
    %%% Response distributions...
    Hist.Y1 = ksdensity(Y1, densityPts, 'kernel', 'Normal');
    Hist.Y2 = ksdensity(Y2, densityPts, 'kernel', 'Normal');
    Hist.Y1c = ksdensity(Y1c, densityPts, 'kernel', 'Normal');
    Hist.Y1m = ksdensity(Y1m, densityPts, 'kernel', 'Normal');
    Hist.Y2m = ksdensity(Y2m, densityPts, 'kernel', 'Normal');
    Hist.Y1b = ksdensity(Y1b, densityPts, 'kernel', 'Normal');
    
    %%% Measure performance and display...
    clearvars ERR
    ERR.Matrix(1, :) = [NRMSE(Y1, Y1m), NMAE(Y1, Y1m), corr(Y1, Y1m, 'type', 'spearman')];
    ERR.Matrix(2, :) = [NRMSE(Y1, Y2m), NMAE(Y1, Y2m), corr(Y1, Y2m, 'type', 'spearman')];
    ERR.Matrix(3, :) = [NRMSE(Y1, Y1c), NMAE(Y1, Y1c), corr(Y1, Y1c, 'type', 'spearman')];
    ERR.Matrix(4, :) = [NRMSE(Y1, Y1b), NMAE(Y1, Y1b), corr(Y1, Y1b, 'type', 'spearman')];
    ERR.Table = array2table(round(ERR.Matrix', 4), 'variablenames', PredApp, 'rownames', {'NRMSE', 'NMAE', 'SCC'});
%     fprintf(1, 'Prediction performance = \n'),            disp(ERR.Table)
    
    % Prediction histograms...
    figure(34),           clf
    subplot(231),        histogram(Y1, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1, 'r-.', 'linewidth', 2),      title(sprintf('Target %s Response', dsFlag{1}))
    subplot(234),       histogram(Y2, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y2, 'r-.', 'linewidth', 2),     title(sprintf('%s Response', dsFlag{2}))
    subplot(235),       histogram(Y1b, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1b, 'r-.', 'linewidth', 2),    title(sprintf('%s Model Prediction', dsFlag{2}))
    subplot(232),       histogram(Y1m, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1m, 'r-.', 'linewidth', 2),    title('{\color[rgb]{0,0.3,0.8}DMTL} Prediction')
    subplot(233),       histogram(Y1c, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1c, 'r-.', 'linewidth', 2),    title('{\color[rgb]{0,0.3,0.8}CORAL-TL} Prediction')
    suptitle([{sprintf('\\bfPrediction of Tumor Culture Response from Cell line Data, Biomarker = {\\color[rgb]{0,0.5,0.4}%s}',...
                        biomarkers{chosenBMidx})}; {'{\color[rgb]{0.7,0,0}Red Line}: Kernel Density Estimate'}])
    % % saveas(figure(34), ['./Figures__MB_CCLE_GDSC_Analysis/',...
    % %                                                         'Prediction_TLP_CP_MCG_', biomarkers{chosenBMidx}, '.fig'])
    % % saveas(figure(34), ['./Figures__MB_CCLE_GDSC_Analysis/',...
    % %                                                         'Prediction_TLP_CP_MCG_', biomarkers{chosenBMidx}, '.png'])
    % saveas(figure(34), [SAVEPATH, sprintf('TLP_SP_Prediction_%s.png', biomarkers{chosenBMidx})])
    % saveas(figure(34), [SAVEPATH, sprintf('TLP_SP_Prediction_MC_%s.png', biomarkers{chosenBMidx})])
    
    % Prediction distributions...
    % Histograms...
    figure(35),           clf
    subplot(231),        histogram(Y1, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1, 'r-.', 'linewidth', 2),      title(sprintf('Target %s Response', dsFlag{1}))
    subplot(234),       histogram(Y2, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y2, 'r-.', 'linewidth', 2),      title(sprintf('%s Response', dsFlag{2}))
    subplot(236),       histogram(Y1b, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1b, 'r-.', 'linewidth', 2),     title(sprintf('%s Model Prediction', dsFlag{2}))
    subplot(232),       histogram(Y2m, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y2m, 'r-.', 'linewidth', 2),    title('TL_{SS} Prediction')
    subplot(233),       histogram(Y1m, numBin, 'normalization', 'pdf')
    hold on,                plot(densityPts, Hist.Y1m, 'r-.', 'linewidth', 2),     title('TL Prediction')
    suptitle([{sprintf('\\bfComparison of Predictions for Tumor data, Biomarker = {\\color[rgb]{0,0.5,0.4}%s}',...
                        biomarkers{chosenBMidx})}; {'{\color[rgb]{0.5,0,0}Red Line}: Kernel density estimate'}])
    % saveas(figure(35), [SAVEPATH, sprintf('TLP_SP_Prediction2_%s.png', biomarkers{chosenBMidx})])
    % saveas(figure(35), [SAVEPATH, sprintf('TLP_SP_Prediction2_MC_%s.png', biomarkers{chosenBMidx})])
    
    %%% Complete Tables...
    RESULTS.NRMSE{chosenBMidx, :} = ERR.Table{1, :};        RESULTS.NMAE{chosenBMidx, :} = ERR.Table{2, :};
    RESULTS.SCC{chosenBMidx, :} = ERR.Table{3, :};             RESULTS.numFeat(chosenBMidx) = m;
    
%     fprintf('\n')
%     disp(array2table(round([RESULTS.NRMSE{chosenBMidx, :}; RESULTS.NMAE{chosenBMidx, :};... 
%                                             RESULTS.SCC{chosenBMidx, :}], 4), 'RowNames', {'NRMSE', 'NMAE', 'SCC'},... 
%                                             'VariableNames', PredApp))
end
RESULTS.NRMSE{end, :} = mean(RESULTS.NRMSE{1:q, :}, 1);
RESULTS.NMAE{end, :} = mean(RESULTS.NMAE{1:q, :}, 1);
RESULTS.SCC{end, :} = mean(RESULTS.SCC{1:q, :}, 1);
RESULTS.numFeat(end) = mean(RESULTS.numFeat(1:q));

fprintf(2, 'Mean performance over %d biomarkers...\n', q)
fprintf(1, '\tm_opt = %d, m_avg = %d\n', m_opt, round(RESULTS.numFeat(end)))
RESULTS.summary = array2table(round([RESULTS.NRMSE{end, :}; RESULTS.NMAE{end, :}; RESULTS.SCC{end, :}], 4),...
                                                            'RowNames', {'NRMSE', 'NMAE', 'SCC'}, 'VariableNames', PredApp);
disp(RESULTS.summary)

