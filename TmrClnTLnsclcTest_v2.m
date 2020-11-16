clc;    clear;      %close all

DATAPATH = sprintf('C:\\Users\\%s\\', getenv('username'));
DATAPATH = {[DATAPATH, 'Dropbox (Texas Tech)\#DataStorageSRD\TCGA_LUAD_cBioPortal\'];
                            [DATAPATH, 'Dropbox (Texas Tech)\#DataStorageSRD\TCGA_LUSC_cBioPortal\'];
                            [DATAPATH, 'Dropbox (Texas Tech)\#DataStorageSRD\CCLE_cBioPortal\cellline_ccle_broad\'];
                            [DATAPATH, 'Google Drive\Data CCLE and GDSC\GDSC\Version 7\']};
FILENAME = {'data_RNA_Seq_v2_expression_median_tcga_luad.txt'; 
                        'data_RNA_Seq_v2_expression_median_tcga_lusc.txt';
                        'data_expression_median_ccle.txt'; 'gdsc_v7_gene_expression_processed_SRD_09_27_2018.xlsx';
                        'data_clinical_sample_ccle.txt'; 'gdsc_v7_Cell_Lines_Details.xlsx'};
EmptyDataFlag = {'.', 'NA', 'N/A'};

%% Data Extraction...
ThZero = 60;

%%% TCGA LUAD (provisional)...
dataTableGX1 = readtable([DATAPATH{1}, FILENAME{1}], 'TreatAsEmpty', EmptyDataFlag);
GX1 = dataTableGX1{:, 3:end}';                            genes1 = dataTableGX1(:, 1:2);
samples1 = dataTableGX1.Properties.VariableNames(3:end)';
nZeros1 = [sum(GX1 == 0, 1); mean(GX1 == 0, 1)*100]';       idxZeros1 = (nZeros1(:, 2) > ThZero);       % Checking zero%
GX1 = GX1(:, ~idxZeros1);                                   genes1 = genes1(~idxZeros1, :);
fprintf(1, '# of LUAD genes discarded due to high zero%% = '),     fprintf(2, '%d\n', sum(idxZeros1))
[~, ind11] = sort(lower(genes1.Hugo_Symbol));     [~, ind12] = sort(lower(samples1));           % Alphabetical sorting
genes1 = genes1(ind11, :);                                     GX1 = GX1(:, ind11);
samples1 = samples1(ind12);                                  GX1 = GX1(ind12, :);

%%% TCGA LUSC (provisional)...
dataTableGX2 = readtable([DATAPATH{2}, FILENAME{2}], 'TreatAsEmpty', EmptyDataFlag);
GX2 = dataTableGX2{:, 3:end}';                           genes2 = dataTableGX2(:, 1:2);
samples2 = dataTableGX2.Properties.VariableNames(3:end)';
nZeros2 = [sum(GX2 == 0, 1); mean(GX2 == 0, 1)*100]';       idxZeros2 = (nZeros2(:, 2) > ThZero);       % Checking zero%
GX2 = GX2(:, ~idxZeros2);                                  genes2 = genes2(~idxZeros2, :);
fprintf(1, '# of LUSC genes discarded due to high zero%% = '),     fprintf(2, '%d\n', sum(idxZeros2))
[~, ind21] = sort(lower(genes2.Hugo_Symbol));    [~, ind22] = sort(lower(samples2));           % Alphabetical sorting
genes2 = genes2(ind21, :);                                    GX2 = GX2(:, ind21);
samples2 = samples2(ind22);                                GX2 = GX2(ind22, :);

% [~, ig1, ig2] = intersect(lower(genes1.Hugo_Symbol), lower(genes2.Hugo_Symbol), 'sorted');
% genes12 = genes1(ig1, :);                                       % Common NSCLC genes

%%% CCLE (vNew)...
dataTableGX3 = readtable([DATAPATH{3}, FILENAME{3}], 'TreatAsEmpty', EmptyDataFlag);
dataGX3 = dataTableGX3{:, 3:end}';                    genes3 = dataTableGX3(:, 1:2);
CLdes3 = dataTableGX3.Properties.VariableDescriptions(3:end)';
idx3 = ~cellfun(@isempty, CLdes3);                     CLdes3 = split(CLdes3(idx3), {' ', ''''});
CL3 = dataTableGX3.Properties.VariableNames(3:end)';           CL3(idx3) = CLdes3(:, end-1);                 % Proper CL names
CLtable3 = readtable([DATAPATH{3}, FILENAME{5}], 'HeaderLines', 4, 'ReadVariableNames', true);   % CL info
[~, ~, idx3] = intersect(CL3, CLtable3.SAMPLE_ID, 'stable');   CLinfo3 = CLtable3(idx3, 2:5);               % Relevant info
tumorType3 = sort(unique(CLinfo3.TUMOR_TYPE));                 lungInd3 = strcmpi(CLinfo3.TUMOR_TYPE, 'lung_nsc');
lungCL3 = CL3(lungInd3);        GX3 = dataGX3(lungInd3, :);     lungCLinfo3 = CLinfo3(lungInd3, :);        % NSCLC data
[~, ind31] = sort(lower(genes3.Hugo_Symbol));    [~, ind32] = sort(lower(lungCL3));           % Alphabetical sorting
genes3 = genes3(ind31, :);                                    GX3 = GX3(:, ind31);
lungCL3 = lungCL3(ind32);                                   GX3 = GX3(ind32, :);                lungCLinfo3 = lungCLinfo3(ind32, :);

%%% GDSC (v7)...
dataTableGX4 = readtable([DATAPATH{4}, FILENAME{4}], 'TreatAsEmpty', EmptyDataFlag);
dataGX4 = dataTableGX4{:, 2:end}';                    genes4 = dataTableGX4.Gene_CL;
CLdes4 = dataTableGX4.Properties.VariableDescriptions(2:end)';
idx4 = ~cellfun(@isempty, CLdes4);                     CLdes4 = split(CLdes4(idx4), {' ', ''''});
CL4 = dataTableGX4.Properties.VariableNames(2:end)';           CL4(idx4) = CLdes4(:, end-1);                 % Proper CL names
CLtable4 = readtable([DATAPATH{4}, FILENAME{6}], 'ReadVariableNames', true);                               % CL info
CLinfo4 = CLtable4(1:end-1, [1, 8:10]);                 CLinfo4.Sample_Name = strrep(CLinfo4.Sample_Name, '-', '');
[~, ~, idx4] = intersect(CL4, join(CLinfo4{:, 1:2}, '_'), 'stable');                                                             % Relevant info
CLinfo4 = CLinfo4(idx4, :);                                  tumorType4 = sort(unique(CLinfo4.GDSC_Tissue_descriptor_1));
lungInd4 = strcmpi(CLinfo4.GDSC_Tissue_descriptor_1, 'lung_nsclc');
lungCL4 = CL4(lungInd4);        GX4 = dataGX4(lungInd4, :);     lungCLinfo4 = CLinfo4(lungInd4, :);        % NSCLC data
[~, ind41] = sort(lower(genes4));                          [~, ind42] = sort(lower(lungCL4));           % Alphabetical sorting
genes4 = genes4(ind41, :);                                    GX4 = GX4(:, ind41);
lungCL4 = lungCL4(ind42);                                   GX4 = GX4(ind42, :);                lungCLinfo4 = lungCLinfo4(ind42, :);

% [~, ig3, ig4] = intersect(lower(genes3.Hugo_Symbol), lower(genes4), 'sorted');
% genes34 = genes3(ig3, :);                                       % Common NSCLC genes


%% Matched data...
matGN = [ ];        [~, matGN(:, 1), matGN(:, 2), matGN(:, 3), matGN(:, 4)] = intersectn(lower(genes1.Hugo_Symbol),... 
                                                        lower(genes2.Hugo_Symbol), lower(genes3.Hugo_Symbol), lower(genes4), 'stable');

% samples1, samples2, lungCL3, lungCL4
geneSet = genes1(matGN(:, 1), :);                          % Common genes betn TCGA, CCLE & GDSC
GXset1 = GX1(:, matGN(:, 1));                               GXset2 = GX2(:, matGN(:, 2));
GXset3 = GX3(:, matGN(:, 3));                              GXset4 = GX4(:, matGN(:, 4));

% % Standardization (z-score)...
% GXset1 = zscore(GXset1, [ ], 1);                            GXset2 = zscore(GXset2, [ ], 1);
% GXset3 = zscore(GXset3, [ ], 1);                           GXset4 = zscore(GXset4, [ ], 1);

% Finding biomarkers in data...
biomarkerLungList = readtable('Biomarker_list_nsclc.csv');
[~, ind] = intersect(lower(geneSet.Hugo_Symbol), lower(biomarkerLungList.Symbol), 'stable');
bmLungIdx = false(numel(geneSet.Hugo_Symbol), 1);       bmLungIdx(ind) = true;
biomarkers = geneSet.Hugo_Symbol(bmLungIdx);             XgeneSet = geneSet(~bmLungIdx, :);

% Biomarker expression as response...
Xdata1 = GXset1(:, ~bmLungIdx);                         Ydata1 = GXset1(:, bmLungIdx);
Xdata2 = GXset2(:, ~bmLungIdx);                        Ydata2 = GXset2(:, bmLungIdx);
Xdata3 = GXset3(:, ~bmLungIdx);                        Ydata3 = GXset3(:, bmLungIdx);
Xdata4 = GXset4(:, ~bmLungIdx);                        Ydata4 = GXset4(:, bmLungIdx);
p = sum(~bmLungIdx);                                          q = sum(bmLungIdx);
n1 = numel(samples1);                                           n2 = numel(samples2);
n3 = numel(lungCL3);                                            n4 = numel(lungCL4);

% Standardization (z-score)...
GXset11 = zscore(GXset1, [ ], 1);                          GXset22 = zscore(GXset2, [ ], 1);
GXset33 = zscore(GXset3, [ ], 1);                         GXset44 = zscore(GXset4, [ ], 1);
Xdata11 = zscore(Xdata1, [ ], 1);                          Ydata11 = zscore(Ydata1, [ ], 1);
Xdata22 = zscore(Xdata2, [ ], 1);                         Ydata22 = zscore(Ydata2, [ ], 1);
Xdata33 = zscore(Xdata3, [ ], 1);                         Ydata33 = zscore(Ydata3, [ ], 1);
Xdata44 = zscore(Xdata4, [ ], 1);                         Ydata44 = zscore(Ydata4, [ ], 1);

%%% Distributions...
numBin = 50;

% Biomarkers...
hb = figure(2);            nPts = 1e3;            clf(figure(hb));          Q = q + mod(q, 2);
mbox = [1:Q/2, 2*Q+1:5*Q/2; Q/2+1:Q, 5*Q/2+1:3*Q; Q+1:3*Q/2, 3*Q+1:7*Q/2; 3*Q/2+1:2*Q, 7*Q/2+1:4*Q];
for k = 1 : q
    [Ypdf1, pts1] = ksdensity(Ydata11(:, k), 'NumPoints', nPts);
    [Ypdf2, pts2] = ksdensity(Ydata22(:, k), 'NumPoints', nPts);
    [Ypdf3, pts3] = ksdensity(Ydata33(:, k), 'NumPoints', nPts);
    [Ypdf4, pts4] = ksdensity(Ydata44(:, k), 'NumPoints', nPts);
    
    figure(hb)                  % NSCLC
    subplot(8,Q/2,mbox(1,k)),     histogram(Ydata11(:, k), numBin, 'normalization', 'pdf')
    hold on,                                 plot(pts1, Ypdf1, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0.5,0,0}TCGA-LUAD'])
    subplot(8,Q/2,mbox(2,k)),     histogram(Ydata22(:, k), numBin, 'normalization', 'pdf')
    hold on,                                 plot(pts2, Ypdf2, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0.5,0.2,0.2}TCGA-LUSC'])
    subplot(8,Q/2,mbox(3,k)),     histogram(Ydata33(:, k), numBin, 'normalization', 'pdf')
    hold on,                                 plot(pts3, Ypdf3, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0.5,0}CCLE-NSCLC'])
    subplot(8,Q/2,mbox(4,k)),     histogram(Ydata44(:, k), numBin, 'normalization', 'pdf')
    hold on,                                 plot(pts4, Ypdf4, 'r-.', 'linewidth', 2),       hold off
    title([biomarkers{k}, ' - \color[rgb]{0,0,0.5}GDSC-NSCLC'])
end
defn = ['Row 1, 4 - {\color[rgb]{0,0,0.8}TCGA-LUAD}, Row 2, 5 - {\color[rgb]{0,0,0.8}TCGA-LUSC}, ',... 
                    'Row 3, 6 - {\color[rgb]{0,0,0.8}CCLE-NSCLC}, Row 4, 8 - {\color[rgb]{0,0,0.8}GDSC-NSCLC}'];
figure(hb),          suptitle([{'\bfBiomarker Distributions for {\color[rgb]{0.5,0,0}NSCLC} CLs'}; {defn}])


%% Pairwise dependency...
% Gene-pair combinations...
bmPairs = nchoosek(1:q, 2);       nPairsBM = size(bmPairs, 1);
bmPairList = strrep(cellstr([num2str(bmPairs(:, 1)), repmat('-', [nPairsBM, 1]), num2str(bmPairs(:, 2))]), ' ', '');

%%% Pairwise Dependency - Correlation...
% Biomarkers...
ycorrmat1 = corr(Ydata1, 'type', 'spearman');                  ycorr1 = 1 - squareform(1 - ycorrmat1, 'tovector')';
ycorrmat2 = corr(Ydata2, 'type', 'spearman');                 ycorr2 = 1 - squareform(1 - ycorrmat2, 'tovector')';
ycorrmat3 = corr(Ydata3, 'type', 'spearman');                 ycorr3 = 1 - squareform(1 - ycorrmat3, 'tovector')';
ycorrmat4 = corr(Ydata4, 'type', 'spearman');                 ycorr4 = 1 - squareform(1 - ycorrmat4, 'tovector')';
ycorr_struct = [corr(ycorr1, ycorr2, 'type', 'pearson'), corr(ycorr1, ycorr3, 'type', 'pearson'),...
                            corr(ycorr1, ycorr4, 'type', 'pearson'), corr(ycorr2, ycorr3, 'type', 'pearson'),...
                            corr(ycorr2, ycorr4, 'type', 'pearson'), corr(ycorr3, ycorr4, 'type', 'pearson')]';
ycorr_struct = [ycorr_struct, round(ycorr_struct, 2)];

% Biomarkers...
figure(3),                        clf
subplot(131),
plot(1:nPairsBM, ycorr1, 'rs--', 1:nPairsBM, ycorr2, 'rd--', 1:nPairsBM, ycorr3, 'go--', 1:nPairsBM, ycorr4, 'bo--')
hold on,                           line([0, nPairsBM+1], [0.5, 0.5], 'color', 'k', 'linestyle', '--'),         xlim([0, nPairsBM+1])
xticksNoted = floor(linspace(1, nPairsBM, 12));         xticks(xticksNoted),        xticklabels(bmPairList(xticksNoted))
xlabel('\bfGene-pairs'),  ylabel('\bfCorrelation'),    legend({'\bfTCGA-LUAD', '\bfTCGA-LUSC', '\bfCCLE', '\bfGDSC'})
title([{'Correlation vector plot: {\color[rgb]{0,0,0.8}NSCLC}'}; 
            {['({\it\color[rgb]{0.5,0,0}r_{TT} = ', num2str(ycorr_struct(1,2)), ', r_{TC} = ', num2str(ycorr_struct(2,2)),... 
                ', r_{TG} = ', num2str(ycorr_struct(3,2)), ', r_{CG} = ', num2str(ycorr_struct(6,2)), '})']}])
subplot(232),                  scatter(ycorr1, ycorr2, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfTCGA-LUAD'),    ylabel('\bfTCGA-LUSC')
title('TCGA-LUAD vs. TCGA-LUSC')
subplot(233),                  scatter(ycorr1, ycorr3, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfTCGA-LUAD'),    ylabel('\bfCCLE-NSCLC')
title('TCGA-LUAD vs. CCLE-NSCLC')
subplot(235),                  scatter(ycorr1, ycorr4, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfTCGA-LUAD'),    ylabel('\bfGDSC-NSCLC')
title('TCGA-LUAD vs. GDSC-NSCLC')
subplot(236),                  scatter(ycorr3, ycorr4, 24, 'bo', 'filled'),                     box on
xticks(-1:0.2:1),               yticks(-1:0.2:1),                 xlabel('\bfCCLE'),                ylabel('\bfGDSC')
title('CCLE vs. GDSC [{\color[rgb]{0,0,0.5}NSCLC}]')
suptitle('\bfComparison of Gene-pair Correlations for {\color[rgb]{0.5,0,0}NSCLC} Biomarkers')


%%%

%% Additional analyses...
norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));

nPts = 1e3;                                   pts = linspace(0, 1, nPts)';
bmPairs = nchoosek(1:q, 2);          nPairsBM = size(bmPairs, 1);         bmPairsNameList = join(biomarkers(bmPairs), ' - ');
Zstates = {'s00', 's01', 's10', 's11'};                    % 4-state discretization
Zpmf1 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Zpmf2 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Zpmf3 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Zpmf4 = array2table(zeros(nPairsBM, 4), 'VariableNames', Zstates);
Th = [0.25, 0.75];                                                % Distribution tail quantiles
asc_met = 'NMI';                                                % Association metric
YpairLo_asc = NaN(nPairsBM, 4);         YpairUp_asc = NaN(nPairsBM, 4);
for b = 1 : nPairsBM
    Ypair1 = norm01(Ydata1(:, bmPairs(b, :)));        Ypair_pdf1 = [ksdensity(Ypair1(:, 1), pts), ksdensity(Ypair1(:, 2), pts)];
    Ypair2 = norm01(Ydata2(:, bmPairs(b, :)));       Ypair_pdf2 = [ksdensity(Ypair2(:, 1), pts), ksdensity(Ypair2(:, 2), pts)];
    Ypair3 = norm01(Ydata3(:, bmPairs(b, :)));       Ypair_pdf3 = [ksdensity(Ypair3(:, 1), pts), ksdensity(Ypair3(:, 2), pts)];
    Ypair4 = norm01(Ydata4(:, bmPairs(b, :)));       Ypair_pdf4 = [ksdensity(Ypair4(:, 1), pts), ksdensity(Ypair4(:, 2), pts)];
    
    % Otsu's thresholding...
    ot1 = [otsuthresh(Ypair_pdf1(:, 1)), otsuthresh(Ypair_pdf1(:, 2))];
    ot2 = [otsuthresh(Ypair_pdf2(:, 1)), otsuthresh(Ypair_pdf2(:, 2))];
    ot3 = [otsuthresh(Ypair_pdf3(:, 1)), otsuthresh(Ypair_pdf3(:, 2))];
    ot4 = [otsuthresh(Ypair_pdf4(:, 1)), otsuthresh(Ypair_pdf4(:, 2))];
    
    % Discretization...
    z1 = ~(Ypair1 <= ot1);        Zpmf1{b, :} = [mean(sum(z1 == [0, 0], 2) == 2), mean(sum(z1 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z1 == [1, 0], 2) == 2), mean(sum(z1 == [1, 1], 2) == 2)];
    z2 = ~(Ypair2 <= ot2);      Zpmf2{b, :} = [mean(sum(z2 == [0, 0], 2) == 2), mean(sum(z2 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z2 == [1, 0], 2) == 2), mean(sum(z2 == [1, 1], 2) == 2)];
    z3 = ~(Ypair3 <= ot3);      Zpmf3{b, :} = [mean(sum(z3 == [0, 0], 2) == 2), mean(sum(z3 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z3 == [1, 0], 2) == 2), mean(sum(z3 == [1, 1], 2) == 2)];
    z4 = ~(Ypair4 <= ot4);      Zpmf4{b, :} = [mean(sum(z4 == [0, 0], 2) == 2), mean(sum(z4 == [0, 1], 2) == 2),... 
                                                                        mean(sum(z4 == [1, 0], 2) == 2), mean(sum(z4 == [1, 1], 2) == 2)];
    
    % Distribution tail threshold quantiles...
    qt1 = quantile(Ypair1, Th);         qt2 = quantile(Ypair2, Th);         qt3 = quantile(Ypair3, Th);         qt4 = quantile(Ypair4, Th);
    
    % Lower & upper dist. tails...
    YpairLo1 = Ypair1((sum(Ypair1 < qt1(1, :), 2) == 2), :);          YpairHi1 = Ypair1((sum(Ypair1 > qt1(2, :), 2) == 2), :);
    YpairLo2 = Ypair2((sum(Ypair2 < qt2(1, :), 2) == 2), :);        YpairHi2 = Ypair2((sum(Ypair2 > qt2(2, :), 2) == 2), :);
    YpairLo3 = Ypair3((sum(Ypair3 < qt3(1, :), 2) == 2), :);        YpairHi3 = Ypair3((sum(Ypair3 > qt3(2, :), 2) == 2), :);
    YpairLo4 = Ypair4((sum(Ypair4 < qt4(1, :), 2) == 2), :);        YpairHi4 = Ypair4((sum(Ypair4 > qt4(2, :), 2) == 2), :);
    
    % Association betn. pairwise distribution tails...
    switch asc_met
        case 'corr'
            % Correlation betn. tails...
            if numel(YpairLo1) > 2,     YpairLo_asc(b, 1) = corr(YpairLo1(:, 1), YpairLo1(:, 2), 'type', 'spearman');       end
            if numel(YpairLo2) > 2,     YpairLo_asc(b, 2) = corr(YpairLo2(:, 1), YpairLo2(:, 2), 'type', 'spearman');     end
            if numel(YpairLo3) > 2,     YpairLo_asc(b, 3) = corr(YpairLo3(:, 1), YpairLo3(:, 2), 'type', 'spearman');     end
            if numel(YpairLo4) > 2,     YpairLo_asc(b, 4) = corr(YpairLo4(:, 1), YpairLo4(:, 2), 'type', 'spearman');     end
            if numel(YpairHi1) > 2,      YpairUp_asc(b, 1) = corr(YpairHi1(:, 1), YpairHi1(:, 2), 'type', 'spearman');     end
            if numel(YpairHi2) > 2,     YpairUp_asc(b, 2) = corr(YpairHi2(:, 1), YpairHi2(:, 2), 'type', 'spearman');    end
            if numel(YpairHi3) > 2,     YpairUp_asc(b, 3) = corr(YpairHi3(:, 1), YpairHi3(:, 2), 'type', 'spearman');    end
            if numel(YpairHi4) > 2,     YpairUp_asc(b, 4) = corr(YpairHi4(:, 1), YpairHi4(:, 2), 'type', 'spearman');    end
            defn = 'Correlation';
        case 'NMI'
            % Normalized Mutual Info betn. tails...
            %%% [MI, Idx] = MutualInformation(X, opt, scale, dtype, C)
            if numel(YpairLo1) > 2,     YpairLo_asc(b, 1) = MutualInformation(YpairLo1, 'single', 'corr');       end
            if numel(YpairLo2) > 2,     YpairLo_asc(b, 2) = MutualInformation(YpairLo2, 'single', 'corr');      end
            if numel(YpairLo3) > 2,     YpairLo_asc(b, 3) = MutualInformation(YpairLo3, 'single', 'corr');      end
            if numel(YpairLo4) > 2,     YpairLo_asc(b, 4) = MutualInformation(YpairLo4, 'single', 'corr');      end
            if numel(YpairHi1) > 2,      YpairUp_asc(b, 1) = MutualInformation(YpairHi1, 'single', 'corr');      end
            if numel(YpairHi2) > 2,     YpairUp_asc(b, 2) = MutualInformation(YpairHi2, 'single', 'corr');      end
            if numel(YpairHi3) > 2,     YpairUp_asc(b, 3) = MutualInformation(YpairHi3, 'single', 'corr');      end
            if numel(YpairHi4) > 2,     YpairUp_asc(b, 4) = MutualInformation(YpairHi4, 'single', 'corr');      end
            defn = 'Normalized MI';
    end
    
    % Printing results...
    if ~mod(b, 10)
        pmfDiff = array2table(round([Zpmf1{b, :} - Zpmf2{b, :}; Zpmf1{b, :} - Zpmf3{b, :}; Zpmf1{b, :} - Zpmf4{b, :};... 
                                                       Zpmf2{b, :} - Zpmf3{b, :}; Zpmf2{b, :} - Zpmf4{b, :}; Zpmf3{b, :} - Zpmf4{b, :}], 4),...
                                                    'RowNames', {'LA-LS', 'LA-C', 'LA-G', 'LS-C', 'LS-G', 'C-G'}, 'VariableNames', Zstates);
        fprintf(1, '\nBiomarker pair [%d] = ', b),        fprintf(2, '%s\n', bmPairsNameList{b})
        fprintf(1, 'Error betn discrete PMFs = \n'),	disp(pmfDiff)
        ascTail = array2table(round([YpairLo_asc(b, :); YpairUp_asc(b, :)], 4), 'RowNames', {'Lower tail'; 'Upper tail'},...
                                                        'VariableNames', {'TCGA_LUAD', 'TCGA_LUSC', 'CCLE', 'GDSC'});
        fprintf(1, '%s between lower & upper distribution tails = \n', asc_met),        disp(ascTail)
    end
end

%%%
% Visualization...
numBin = 50;
bmPairList = strrep(cellstr([num2str(bmPairs(:, 1)), repmat('-', [nPairsBM, 1]), num2str(bmPairs(:, 2))]), ' ', '');

figure(4),         clf
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
xlabel(['\bf', biomarkers{bmPairs(b, 1)}]),      ylabel(['\bf', biomarkers{bmPairs(b, 2)}])
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}TCGA-LUAD}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')

figure(5),         clf
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
xlabel(['\bf', biomarkers{bmPairs(b, 1)}]),      ylabel(['\bf', biomarkers{bmPairs(b, 2)}])
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}TCGA-LUSC}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')

figure(6),         clf
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
xlabel(['\bf', biomarkers{bmPairs(b, 1)}]),      ylabel(['\bf', biomarkers{bmPairs(b, 2)}])
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}CCLE}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')

figure(7),         clf
ax1 = axes('Position', [0.05, 0.14, 0.15, 0.6]);       ax2 = axes('Position', [0.21, 0.75, 0.65, 0.2]);
ax3 = axes('Position', [0.21, 0.14, 0.65, 0.6], 'YAxisLocation', 'right');
histogram(ax1, Ypair4(:, 2), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax1, [-90, 90]),	hold on
histogram(ax2, Ypair4(:, 1), numBin, 'Normalization', 'pdf', 'FaceColor', [0.8, 0.5, 0]),	view(ax2, [0, 90]),      hold on
scatter(ax3, Ypair4(:, 1), Ypair4(:, 2), 18, 'filled', 'o'),       box on
hold on,       line(ax3, [0, 1], [ot4(2), ot4(2)], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt4(1, 2), qt4(1, 2)], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt4(1, 1), qt4(1, 1)], [0, 1], 'color', [0.4, 0.6, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [ot4(1), ot4(1)], [0, 1], 'color', 'r', 'linewidth', 2, 'linestyle', '--')
hold on,       line(ax3, [0, 1], [qt4(2, 2), qt4(2, 2)], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
hold on,       line(ax3, [qt4(2, 1), qt4(2, 1)], [0, 1], 'color', [0.6, 0.4, 0], 'linewidth', 1.5, 'linestyle', '--')
legend('\bf\itData', '\bfTh_{OTSU}', '\bfTh_{LO}', '\bfTh_{UP}')
xlabel(['\bf', biomarkers{bmPairs(b, 1)}]),      ylabel(['\bf', biomarkers{bmPairs(b, 2)}])
title('\bfPairwise Biomarker Scatter: {\color[rgb]{0,0,0.8}GDSC}', 'Position', [0.4, 1.35, 0], 'FontSize', 14)
set(ax1, 'visible', 'off'),     set(ax2, 'visible', 'off')

%%%
% Linear model fit...
model13 = {fitlm(Zpmf1{:, 1}, Zpmf3{:, 1}), fitlm(Zpmf1{:, 2}, Zpmf3{:, 2}),...
                        fitlm(Zpmf1{:, 3}, Zpmf3{:, 3}), fitlm(Zpmf1{:, 4}, Zpmf3{:, 4})};
model14 = {fitlm(Zpmf1{:, 1}, Zpmf4{:, 1}), fitlm(Zpmf1{:, 2}, Zpmf4{:, 2}),... 
                        fitlm(Zpmf1{:, 3}, Zpmf4{:, 3}), fitlm(Zpmf1{:, 4}, Zpmf4{:, 4})};
model23 = {fitlm(Zpmf2{:, 1}, Zpmf3{:, 1}), fitlm(Zpmf2{:, 2}, Zpmf3{:, 2}),...
                        fitlm(Zpmf2{:, 3}, Zpmf3{:, 3}), fitlm(Zpmf2{:, 4}, Zpmf3{:, 4})};
model24 = {fitlm(Zpmf2{:, 1}, Zpmf4{:, 1}), fitlm(Zpmf2{:, 2}, Zpmf4{:, 2}),... 
                        fitlm(Zpmf2{:, 3}, Zpmf4{:, 3}), fitlm(Zpmf2{:, 4}, Zpmf4{:, 4})};
w13 = [model13{1}.Coefficients.Estimate, model13{2}.Coefficients.Estimate,... 
                model13{3}.Coefficients.Estimate, model13{4}.Coefficients.Estimate];
w14 = [model14{1}.Coefficients.Estimate, model14{2}.Coefficients.Estimate,... 
                model14{3}.Coefficients.Estimate, model14{4}.Coefficients.Estimate];
w23 = [model23{1}.Coefficients.Estimate, model23{2}.Coefficients.Estimate,... 
                model23{3}.Coefficients.Estimate, model23{4}.Coefficients.Estimate];
w24 = [model24{1}.Coefficients.Estimate, model24{2}.Coefficients.Estimate,... 
                model24{3}.Coefficients.Estimate, model24{4}.Coefficients.Estimate];
Wm1 = [mean([w13(1, :); w14(1, :)], 1); mean([w13(2, :); w14(2, :)], 1)];                                % Mean coefficient
xval1 = linspace(0, round(max(max(Zpmf1{:, :})), 1), 100)';       Xval1 = [ones(size(xval1)), xval1];
yfit1 = [Xval1 * Wm1(:, 1), Xval1 * Wm1(:, 2), Xval1 * Wm1(:, 3), Xval1 * Wm1(:, 4)];            % Fitted line
Wm2 = [mean([w23(1, :); w24(1, :)], 1); mean([w23(2, :); w24(2, :)], 1)];                             % Mean coefficient
xval2 = linspace(0, round(max(max(Zpmf2{:, :})), 1), 100)';       Xval2 = [ones(size(xval2)), xval2];
yfit2 = [Xval2 * Wm2(:, 1), Xval2 * Wm2(:, 2), Xval2 * Wm2(:, 3), Xval2 * Wm2(:, 4)];       % Fitted line

% Discretized biomarker-pair PMFs...
figure(8),            clf
subplot(241),       scatter(Zpmf1{:, 1}, Zpmf3{:, 1}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 1}, Zpmf4{:, 1}, 24, 'ro', 'filled'),                    box on
hold on,               plot(xval1, yfit1(:, 1), 'g-'),          xlabel('\bfTumor Culture'),       ylabel('\bfCell line')
legend('\bfLUAD - CCLE', '\bfLUAD - GDSC', 'Location', 'nw'),          title('State = \color[rgb]{0,0,0.8}00')
subplot(242),       scatter(Zpmf1{:, 2}, Zpmf3{:, 2}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 2}, Zpmf4{:, 2}, 24, 'ro', 'filled'),                   box on
hold on,               plot(xval1, yfit1(:, 2), 'g-'),          xlabel('\bfTumor Culture'),      ylabel('\bfCell line')
legend('\bfLUAD - CCLE', '\bfLUAD - GDSC', 'Location', 'nw'),         title('State = \color[rgb]{0,0,0.8}01')
subplot(243),       scatter(Zpmf1{:, 3}, Zpmf3{:, 3}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 3}, Zpmf4{:, 3}, 24, 'ro', 'filled'),                  box on
hold on,               plot(xval1, yfit1(:, 3), 'g-'),          xlabel('\bfTumor Culture'),     ylabel('\bfCell line')
legend('\bfLUAD - CCLE', '\bfLUAD - GDSC', 'Location', 'nw'),        title('State = \color[rgb]{0,0,0.8}10')
subplot(244),       scatter(Zpmf1{:, 4}, Zpmf3{:, 4}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf1{:, 4}, Zpmf4{:, 4}, 24, 'ro', 'filled'),                 box on
hold on,               plot(xval1, yfit1(:, 4), 'g-'),          xlabel('\bfTumor Culture'),    ylabel('\bfCell line')
legend('\bfLUAD - CCLE', '\bfLUAD - GDSC', 'Location', 'nw'),        title('State = \color[rgb]{0,0,0.8}11')
subplot(245),       scatter(Zpmf2{:, 1}, Zpmf3{:, 1}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf2{:, 1}, Zpmf4{:, 1}, 24, 'ro', 'filled'),                    box on
hold on,               plot(xval2, yfit2(:, 1), 'g-'),          xlabel('\bfTumor Culture'),       ylabel('\bfCell line')
legend('\bfLUSC - CCLE', '\bfLUSC - GDSC', 'Location', 'nw'),          title('State = \color[rgb]{0,0,0.8}00')
subplot(246),       scatter(Zpmf2{:, 2}, Zpmf3{:, 2}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf2{:, 2}, Zpmf4{:, 2}, 24, 'ro', 'filled'),                   box on
hold on,               plot(xval2, yfit2(:, 2), 'g-'),          xlabel('\bfTumor Culture'),      ylabel('\bfCell line')
legend('\bfLUSC - CCLE', '\bfLUSC - GDSC', 'Location', 'nw'),         title('State = \color[rgb]{0,0,0.8}01')
subplot(247),       scatter(Zpmf2{:, 3}, Zpmf3{:, 3}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf2{:, 3}, Zpmf4{:, 3}, 24, 'ro', 'filled'),                  box on
hold on,               plot(xval2, yfit2(:, 3), 'g-'),          xlabel('\bfTumor Culture'),     ylabel('\bfCell line')
legend('\bfLUSC - CCLE', '\bfLUSC - GDSC', 'Location', 'nw'),        title('State = \color[rgb]{0,0,0.8}10')
subplot(248),       scatter(Zpmf2{:, 4}, Zpmf3{:, 4}, 24, 'bo', 'filled')
hold on,               scatter(Zpmf2{:, 4}, Zpmf4{:, 4}, 24, 'ro', 'filled'),                 box on
hold on,               plot(xval2, yfit2(:, 4), 'g-'),          xlabel('\bfTumor Culture'),    ylabel('\bfCell line')
legend('\bfLUSC - CCLE', '\bfLUSC - GDSC', 'Location', 'nw'),        title('State = \color[rgb]{0,0,0.8}11')
suptitle(['\bfComparison of {\color[rgb]{0,0,0.8}Discretized} Gene-pair Fractions for Biomarkers ',...
                    '[{\color[rgb]{0.6,0.2,0.2}NSCLC}]'])


% Corr. betn. tails...
figure(9),        clf
subplot(241),   scatter(YpairLo_asc(:, 1), YpairLo_asc(:, 2), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 1))+0.1, -0.1, max(YpairLo_asc(:, 2))+0.1])
xlabel('\bfTCGA-LUAD'),       ylabel('\bfTCGA-LUSC'),       title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(242),   scatter(YpairLo_asc(:, 1), YpairLo_asc(:, 3), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 1))+0.1, -0.1, max(YpairLo_asc(:, 3))+0.1])
xlabel('\bfTCGA-LUAD'),       ylabel('\bfCCLE'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(243),   scatter(YpairLo_asc(:, 1), YpairLo_asc(:, 4), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 1))+0.1, -0.1, max(YpairLo_asc(:, 3))+0.1])
xlabel('\bfTCGA-LUAD'),       ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(244),   scatter(YpairLo_asc(:, 3), YpairLo_asc(:, 4), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairLo_asc(:, 2))+0.1, -0.1, max(YpairLo_asc(:, 3))+0.1])
xlabel('\bfCCLE'),                  ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Lower Tails', defn))
subplot(245),   scatter(YpairUp_asc(:, 1), YpairUp_asc(:, 2), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 1))+0.1, -0.1, max(YpairUp_asc(:, 2))+0.1])
xlabel('\bfTCGA-LUAD'),       ylabel('\bfTCGA-LUSC'),       title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
subplot(246),   scatter(YpairUp_asc(:, 1), YpairUp_asc(:, 3), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 1))+0.1, -0.1, max(YpairUp_asc(:, 3))+0.1])
xlabel('\bfTCGA-LUAD'),       ylabel('\bfCCLE'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
subplot(247),   scatter(YpairUp_asc(:, 1), YpairUp_asc(:, 4), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 1))+0.1, -0.1, max(YpairUp_asc(:, 3))+0.1])
xlabel('\bfTCGA-LUAD'),       ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
subplot(248),   scatter(YpairUp_asc(:, 3), YpairUp_asc(:, 4), 24, 'o', 'filled'),     box on,      %grid on
axis([-0.1, max(YpairUp_asc(:, 2))+0.1, -0.1, max(YpairUp_asc(:, 3))+0.1])
xlabel('\bfCCLE'),                  ylabel('\bfGDSC'),      title(sprintf('{\\color[rgb]{0,0,0.8}%s}: Upper Tails', defn))
suptitle('\bfComparison of Gene-pair Distribution {\itTails} for Biomarkers [{\color[rgb]{0.6,0.2,0.2}NSCLC}]')

% Correlation betn. association vectors...
kLo12 = ~(isnan(YpairLo_asc(:, 1)) | isnan(YpairLo_asc(:, 2)));     kLo13 = ~(isnan(YpairLo_asc(:, 1)) | isnan(YpairLo_asc(:, 3)));
kLo14 = ~(isnan(YpairLo_asc(:, 1)) | isnan(YpairLo_asc(:, 4)));     kLo23 = ~(isnan(YpairLo_asc(:, 2)) | isnan(YpairLo_asc(:, 3)));
kLo24 = ~(isnan(YpairLo_asc(:, 2)) | isnan(YpairLo_asc(:, 4)));    kLo34 = ~(isnan(YpairLo_asc(:, 3)) | isnan(YpairLo_asc(:, 4)));
kUp12 = ~(isnan(YpairUp_asc(:, 1)) | isnan(YpairUp_asc(:, 2)));   kUp13 = ~(isnan(YpairUp_asc(:, 1)) | isnan(YpairUp_asc(:, 3)));
kUp14 = ~(isnan(YpairUp_asc(:, 1)) | isnan(YpairUp_asc(:, 4)));   kUp23 = ~(isnan(YpairUp_asc(:, 2)) | isnan(YpairUp_asc(:, 3)));
kUp24 = ~(isnan(YpairUp_asc(:, 2)) | isnan(YpairUp_asc(:, 4)));  kUp34 = ~(isnan(YpairUp_asc(:, 3)) | isnan(YpairUp_asc(:, 4)));
YpairTail_ascStruct = array2table(round([corr(YpairLo_asc(kLo12, 1), YpairLo_asc(kLo12, 2), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo13, 1), YpairLo_asc(kLo13, 3), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo14, 1), YpairLo_asc(kLo14, 4), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo23, 2), YpairLo_asc(kLo23, 3), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo24, 2), YpairLo_asc(kLo24, 4), 'type', 'spearman'),...
                                            corr(YpairLo_asc(kLo34, 3), YpairLo_asc(kLo34, 4), 'type', 'spearman');...
                                            corr(YpairUp_asc(kUp12, 1), YpairUp_asc(kUp12, 2), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp13, 1), YpairUp_asc(kUp13, 3), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp14, 1), YpairUp_asc(kUp14, 4), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp23, 2), YpairUp_asc(kUp23, 3), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp24, 2), YpairUp_asc(kUp24, 4), 'type', 'spearman'),...
                                            corr(YpairUp_asc(kUp34, 3), YpairUp_asc(kUp34, 4), 'type', 'spearman')], 4),...
                                        'rownames', {'Lower tail'; 'Upper tail'}, 'variablename',... 
                                        {'LA_LS', 'LA_C', 'LA_G', 'LS_C', 'LS_G', 'C_G'});
disp(YpairTail_ascStruct)


%% Feature selection - TIME-CONSUMING LOOP!!!...
% RELIEFF gene ranking...
K = 10;
rank1 = zeros(p, q);        rank2 = zeros(p, q);
rank3 = zeros(p, q);        rank4 = zeros(p, q);
tic
for k = 1 : q
    fprintf('Biomarker %d: %s\n', k, biomarkers{k})
    rank1(:, k) = relieff(Xdata1, Ydata1(:, k), K, 'method', 'regression');
    rank2(:, k) = relieff(Xdata2, Ydata2(:, k), K, 'method', 'regression');
    rank3(:, k) = relieff(Xdata3, Ydata3(:, k), K, 'method', 'regression');
    rank4(:, k) = relieff(Xdata4, Ydata4(:, k), K, 'method', 'regression');
end
FEATFILENAME = 'relieff_ranked_genes_nsclc_TCG.mat';
save(FEATFILENAME, 'rank1', 'rank2', 'rank3', 'rank4')
toc


%% Transfer Learning...
% Load RELEIFF features...
FEATFILENAME = 'relieff_ranked_genes_nsclc_TCG.mat';
load(FEATFILENAME, 'rank1', 'rank2', 'rank3', 'rank4')

% Metric function definitions...
norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));       % Normalize in [0, 1]
NMAE = @(y, y_hat) mean(abs(y - y_hat)) / mean(abs(y - mean(y)));
NRMSE = @(y, y_hat) sqrt(mean((y - y_hat).^2)) / std(y, 1);
RSQ = @(y, y_hat) 1 - mean((y - y_hat).^2) / var(y, 1);

fprintf(2, 'Tumor cultures to cell line data Transfer Learning...\n')
PredApp = {'DMTL', 'DMTL_SS', 'CATL', 'BL'};
RESULTS = array2table(zeros(q+1, numel(PredApp)), 'RowNames', [biomarkers; {'Mean'}], 'VariableNames', PredApp);
RESULTS = struct('NRMSE', RESULTS, 'NMAE', RESULTS, 'SCC', RESULTS, 'numFeat', zeros(q+1, 1));
qn = 1 : q;
for chosenBMidx = qn
    %%% Feature selection...
    fprintf(1, 'Chosen biomarker = '),       fprintf(2, '%s\n', biomarkers{chosenBMidx})
    ranks = [rank1(:, chosenBMidx), rank2(:, chosenBMidx), rank3(:, chosenBMidx), rank4(:, chosenBMidx)];
    
    % Tumor-to-CL TL case...
    % fprintf(1, 'Extracting common CL features...\n')                    % Extracting a common set of size m_opt
    nGN = 300;                  m_opt = 150;        gnRank = intersectn(ranks(1:nGN, 3), ranks(1:nGN, 4), 'sorted');
    m = numel(gnRank);       nI = 0;                 m0 = m;
    while m < m_opt
        nI = nI + 1;               nGN = nGN + 100;
        gnRank = intersectn(ranks(1:nGN, 3), ranks(1:nGN, 4), 'sorted');        m = numel(gnRank);
    end
    % fprintf(1, '\t#genes used = %d, #Iterations = %d\n\tInitial size = %d, Final size = %d\n', nGN, nI, m0, m)
    
    % CL-to-Tumor TL case...
    % fprintf(1, 'Extracting common tumor culture features...\n')                    % Extracting a common set of size m_opt
    nGN = 300;                                m_opt = 150;          gnRank_rev = intersectn(ranks(1:nGN, 1), ranks(1:nGN, 2), 'sorted');
    m_rev = numel(gnRank_rev);       nI = 0;                   m0_rev = m_rev;
    while m_rev < m_opt
        nI = nI + 1;               nGN = nGN + 100;
        gnRank_rev = intersectn(ranks(1:nGN, 1), ranks(1:nGN, 2), 'sorted');        m_rev = numel(gnRank_rev);
    end
    % fprintf(1, '\t#genes used = %d, #Iterations = %d\n\tInitial size = %d, Final size = %d\n', nGN, nI, m0_rev, m_rev)
    
    % Defining primary & secondary sets...
    fprintf(1, 'Tumor cultures to cell line data Transfer Learning...\n')
    dsChoice = 1;
    switch dsChoice
        case 1                                                              % Tumor => primary, CL => secondary
            dsFlag = {'TCGA-NSCLC', 'CCLE + GDSC'};
            gnList = XgeneSet.Hugo_Symbol(gnRank);                                m = numel(gnList);
            X1 = [Xdata1(:, gnRank); Xdata2(:, gnRank)];                            X2 = [Xdata3(:, gnRank); Xdata4(:, gnRank)];
            Y1 = [(Ydata1(:, chosenBMidx)); (Ydata2(:, chosenBMidx))];          
            Y2 = [(Ydata3(:, chosenBMidx)); (Ydata4(:, chosenBMidx))];
        case 2                                                             % CL => primary, Tumor => secondary
            dsFlag = {'CCLE + GDSC', 'TCGA-NSCLC'};
            gnList = XgeneSet.Huwgo_Symbol(gnRank_rev);                          m = numel(gnList);
            X1 = [Xdata3(:, gnRank_rev); Xdata4(:, gnRank_rev)];               X2 = [Xdata1(:, gnRank_rev); Xdata2(:, gnRank_rev)];
            Y1 = [(Ydata3(:, chosenBMidx)); (Ydata4(:, chosenBMidx))];          
            Y2 = [(Ydata1(:, chosenBMidx)); (Ydata2(:, chosenBMidx))];
        case 3                                                             % Tumor => primary, CCLE => secondary
            dsFlag = {'TCGA-NSCLC', 'CCLE'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 3));                          m = numel(gnList);
            X1 = [Xdata1(:, ranks(1:m_opt, 3)); Xdata2(:, ranks(1:m_opt, 3))];      X2 = Xdata3(:, ranks(1:m_opt, 3));
            Y1 = [Ydata1(:, chosenBMidx); Ydata2(:, chosenBMidx)];                   Y2 = Ydata3(:, chosenBMidx);
        case 4                                                             % CCLE => primary, Tumor => secondary
            dsFlag = {'CCLE', 'TCGA-NSCLC'};
            gnList = XgeneSet.Hugo_Symbol(gnRank_rev);          m = numel(gnList);
            X1 = Xdata3(:, gnRank_rev);                                     X2 = [Xdata1(:, gnRank_rev); Xdata2(:, gnRank_rev)];
            Y1 = Ydata3(:, chosenBMidx);                                   Y2 = [Ydata1(:, chosenBMidx); Ydata2(:, chosenBMidx)];
        case 5                                                             % Tumor => primary, GDSC => secondary
            dsFlag = {'TCGA-NSCLC', 'GDSC'};
            gnList = XgeneSet.Hugo_Symbol(ranks(1:m_opt, 4));                          m = numel(gnList);
            X1 = [Xdata1(:, ranks(1:m_opt, 4)); Xdata2(:, ranks(1:m_opt, 4))];      X2 = Xdata4(:, ranks(1:m_opt, 4));
            Y1 = [Ydata1(:, chosenBMidx); Ydata2(:, chosenBMidx)];                   Y2 = Ydata4(:, chosenBMidx);
        case 6                                                             % GDSC => primary, Tumor => secondary
            dsFlag = {'GDSC', 'TCGA-NSCLC'};
            gnList = XgeneSet.Hugo_Symbol(gnRank_rev);          m = numel(gnList);
            X1 = Xdata4(:, gnRank_rev);                                     X2 = [Xdata1(:, gnRank_rev); Xdata2(:, gnRank_rev)];
            Y1 = Ydata4(:, chosenBMidx);                                   Y2 = [Ydata1(:, chosenBMidx); Ydata2(:, chosenBMidx)];
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
    datadir = sprintf('NSCLC_%s_%s', dsChar{3}, dsChar{2});
    filename = sprintf('Data_BM_%d_f%d_%s%s_NSCLC.mat', chosenBMidx, m_opt, dsChar{2}, dsChar{1});
    save(fullfile(datapath, datadir, filename), 'X1', 'X2', 'Y1', 'Y2')
    
    % Distributions...
    figure(10)
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
    
    %%% Transfer Learning predictions...
    [Y1c, X2a, CorrMatCA] = CorrAlignTransLearn(X1, X2, Y2);
    [Y1m, Y2m, Map, densityPts, CorrMatDM] = DistMatchTransLearn(X1, Y1, X2, Y2);
    
    
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
    
    % Prediction histograms...
    figure(11),             clf
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
    
    % Prediction distributions...
    % Histograms...
    figure(12),            clf
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
    
    %%% Complete Tables...
    RESULTS.NRMSE{chosenBMidx, :} = ERR.Table{1, :};        RESULTS.NMAE{chosenBMidx, :} = ERR.Table{2, :};
    RESULTS.SCC{chosenBMidx, :} = ERR.Table{3, :};             RESULTS.numFeat(chosenBMidx) = m;
    
%     fprintf('\n')
%     disp(array2table(round([RESULTS.NRMSE{chosenBMidx, :}; RESULTS.NMAE{chosenBMidx, :};...
%                                                 RESULTS.SCC{chosenBMidx, :}], 4), 'RowNames', {'NRMSE', 'NMAE', 'SCC'},...
%                                                 'VariableNames', PredApp))
end
RESULTS.NRMSE{end, :} = mean(RESULTS.NRMSE{qn, :}, 1);
RESULTS.NMAE{end, :} = mean(RESULTS.NMAE{qn, :}, 1);
RESULTS.SCC{end, :} = mean(RESULTS.SCC{qn, :}, 1);
RESULTS.numFeat(end) = mean(RESULTS.numFeat(qn));

fprintf(2, 'Mean performance over %d biomarkers...\n', length(qn))
fprintf(1, '\tm_opt = %d, m_avg = %d\n', m_opt, round(RESULTS.numFeat(end)))
RESULTS.summary = array2table(round([RESULTS.NRMSE{end, :}; RESULTS.NMAE{end, :}; RESULTS.SCC{end, :}], 4),...
                                                            'RowNames', {'NRMSE', 'NMAE', 'SCC'}, 'VariableNames', PredApp);
disp(RESULTS.summary)











