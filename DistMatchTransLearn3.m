function [Ypred_PS, Ypred_SS, Map, distPts, CorrStr, X2m, XTm, Hist] = DistMatchTransLearn3(Xp, Yp, Xs, Ys, Xt)
    % Initial checks....
    if (size(Xs, 1) ~= length(Ys)),          error('Secondary set must have matched datasets!');     end
    if (size(Xp, 2) ~= size(Xs, 2)),        error('Number of features must be same!');                   end
    if (size(Xp, 2) ~= size(Xt, 2)),        error('Number of features must be same!');                   end
    
    % Assign datasets...
    %norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));       % Normalize in [0, 1]
    %X1 = norm01(Xp);        Y1 = (Yp(:));        n1 = size(X1, 1);
    %X2 = norm01(Xs);        Y2 = (Ys(:));        n2 = size(X2, 1);
    %XT = norm01(Xt);                                   nT = size(XT, 1);
    
    n1 = size(Xp, 1);     n2 = size(Xs, 1);      [nT, m] = size(Xt);
    stat_p = struct('min', min(Xp, [ ], 1), 'max', max(Xp, [ ], 1), 'range', range(Xp, 1), 'mu', mean(Xp, 1), 'sd', std(Xp, 0, 1));
    stat_s = struct('min', min(Xs, [ ], 1), 'max', max(Xs, [ ], 1), 'range', range(Xs, 1), 'mu', mean(Xs, 1), 'sd', std(Xs, 0, 1));
    
    norm_switch = 1;
    switch norm_switch
        case 1
            X1 = (Xp - ones(n1, 1) * stat_p.min) ./ (ones(n1, 1) * stat_p.range);                   Y1 = Yp(:);
            X2 = (Xs - ones(n2, 1) * stat_s.min) ./ (ones(n2, 1) * stat_s.range);                  Y2 = Ys(:);
            XT = (Xt - ones(nT, 1) * stat_p.min) ./ (ones(nT, 1) * stat_p.range);
        case 2
            X1 = (Xp - ones(n1, 1) * stat_p.mu) ./ (ones(n1, 1) * stat_p.sd);                   Y1 = Yp(:);
            X2 = (Xs - ones(n2, 1) * stat_s.mu) ./ (ones(n2, 1) * stat_s.sd);                  Y2 = Ys(:);
            XT = (Xt - ones(nT, 1) * stat_p.mu) ./ (ones(nT, 1) * stat_p.sd);
        case 3
            X1 = Xp;                   Y1 = Yp(:);
            X2 = Xs;                  Y2 = Ys(:);
            XT = Xt;
    end
    
    % nPts = 1e3;         Y_distPts = [linspace(min(Y1), max(Y1), nPts); linspace(min(Y2), max(Y2), nPts)]';
    
    tic
    %%% Transfer Learning w/ Histogram Matching...
    % Initialize parameters...
	nPts = 500;         %distPts = linspace(0, 1, nPts+1)';                       % Sample points
%     distPts = struct('X', linspace(min([X1(:); X2(:)]), max([X1(:); X2(:)]) , nPts+1),...
%                                                             'Y', linspace(min([Y1(:); Y2(:)]), max([Y1(:); Y2(:)]), nPts+1));
    distPts = struct('X', linspace(0, 1, nPts), 'Y', linspace(0, 1, nPts));
    MapCalc = 0;      L = 256 * (MapCalc == 0) + nPts * (MapCalc ~= 0);
    data2idx = @(z, zmax) round((L - 1) * z / zmax) + 1;         % Convert to histogram map indices
    Hist = struct('X1', zeros(nPts, m), 'X2', zeros(nPts, m), 'X2m', zeros(nPts, m), 'X1m', zeros(nPts, m),... 
                            'XT', zeros(nPts, m), 'XTm', zeros(nPts, m),...
                            'Y1', zeros(nPts, 1), 'Y2', zeros(nPts, 1), 'Y2m', zeros(nPts, 1), 'Y1m', zeros(nPts, 1));
    Map = struct('X', zeros(L, m), 'Y', zeros(L, 1), 'Yinv', zeros(L, 1));
    
    %%% Covariate Map [DS1 => DS2]: per feature matching...
    X2m = zeros(n1, m);         XTm = zeros(nT, m);
    for j = 1 : m
        Hist.X2(:, j) = ksdensity(X2(:, j), distPts.X, 'kernel', 'Normal');
        Hist.X1(:, j) = ksdensity(X1(:, j), distPts.X, 'kernel', 'Normal');
        Hist.XT(:, j) = ksdensity(XT(:, j), distPts.X, 'kernel', 'Normal');
                
%         Hist.X2(:, j) = histcounts(X2(:, j), distPts.X, "Normalization", "pdf") / nPts;
%         Hist.X1(:, j) = histcounts(X1(:, j), distPts.X, "Normalization", "pdf") / nPts;
%         Hist.XT(:, j) = histcounts(XT(:, j), distPts.X, "Normalization", "pdf") / nPts;
        
        % Mapping for j-th feature...
        switch MapCalc
            case 0
                [~, Map.X(:, j)] = histeq(X1(:, j), Hist.X2(:, j));
            case 1
                Map.X(:, j) = histmatch(Hist.X1(:, j), distPts, Hist.X2(:, j), distPts);
        end
        X2m(:, j) = Map.X(data2idx(X1(:, j), 1), j);                        % Mapped primary X
        Hist.X2m(:, j) = ksdensity(X2m(:, j), distPts.X, 'kernel', 'Normal');
%         Hist.X2m(:, j) = histcounts(X2m(:, j), distPts.X, "Normalization", "probability");
        
        XTm(:, j) = Map.X(data2idx(XT(:, j), 1), j);                        % Mapped primary X
        Hist.XTm(:, j) = ksdensity(XTm(:, j), distPts.X, 'kernel', 'Normal');
%         Hist.XTm(:, j) = histcounts(XTm(:, j), distPts.X, "Normalization", "probability");
    end
    
    %%% Response Maps...
%     Hist.Y2 = ksdensity(Y2, distPts.Y, 'kernel', 'Normal');
%     Hist.Y1 = ksdensity(Y1, distPts.Y, 'kernel', 'Normal');
    
    Hist.Y2 = histcounts(Y2, distPts.Y, "Normalization", "pdf") / nPts;
    Hist.Y1 = histcounts(Y1, distPts.Y, "Normalization", "pdf") / nPts;
    
    % Forward & Inverse Maps...
    switch MapCalc
        case 0
            [~, Map.Y] = histeq(Y1, Hist.Y2);            Map.Y = Map.Y(:);              % DS1 => DS2
            [~, Map.Yinv] = histeq(Y2, Hist.Y1);        Map.Yinv = Map.Yinv(:);     % DS2 => DS1
        case 1
            Map.Y = histmatch(Hist.Y1, distPts, Hist.Y2, distPts);                      % DS1 => DS2
            Map.Yinv = histmatch(Hist.Y2, distPts, Hist.Y1, distPts);                  % DS2 => DS1
    end
    Y2m = Map.Y(data2idx(Y1, 1));                                                                  % Mapped primary Y
    Y1m = Map.Yinv(data2idx(Y2, 1));                                                              % Inverse mapped secondary Y
%     Hist.Y2m = ksdensity(Y2m, distPts.Y, 'kernel', 'Normal');
%     Hist.Y1m = ksdensity(Y1m, distPts.Y, 'kernel', 'Normal');
    
    Hist.Y2m = histcounts(Y2m, distPts.Y, "Normalization", "pdf") / nPts;
    Hist.Y1m = histcounts(Y1m, distPts.Y, "Normalization", "pdf") / nPts;
    
    %%% TL prediction...
    ModelSelect = 'EN';
    TrainX = zscore(X2);	TrainY = Y2;        TestX = zscore(XTm);
    switch ModelSelect
        case 'EN'
            rng(0);             Kf = 5;             Alpha = 0.005;
            [Beta, FitInfo] = lasso(TrainX, TrainY, 'alpha', Alpha, 'cv', Kf);          % EN
            beta = Beta(:, FitInfo.IndexMinMSE);              beta0 = FitInfo.Intercept(FitInfo.IndexMinMSE);
            PredY = beta0 + TestX * beta;                          YTp = PredY;         YTp(YTp < 0) = 0;       YTp(YTp > 1) = 1;
        case 'SVR'
            Kernel = 'Gaussian';        PolyOrd = [ ];
            SVR = fitrsvm(TrainX, TrainY, 'KernelFunction', Kernel, 'PolynomialOrder', PolyOrd, 'BoxConstraint', 10);
            PredY = predict(SVR, TestX);         YTp = PredY;        YTp(YTp < 0) = 0;       YTp(YTp > 1) = 1;            
        case 'RF'
            rng(0);             nTree = 200;        RF = TreeBagger(nTree, TrainX, TrainY, 'method', 'regression', 'MinLeafSize', 7);
            PredY = predict(RF, TestX);         YTp = PredY;        YTp(YTp < 0) = 0;       YTp(YTp > 1) = 1;
    end
    Yt_p = Map.Yinv(data2idx(YTp, 1));                                        % Map back to primary space
    %%% TL model finished.
    toc
    
    %%% Outputs...
    Ypred_PS = Yt_p;      Ypred_SS = YTp;
    
    clearvars CorrStr
    CorrStr.Xp = corr(X1, 'type', 'pearson');          CorrStr.Xp1D = 1 - squareform(1 - CorrStr.Xp, 'tovector')';
    CorrStr.Xs = corr(X2, 'type', 'pearson');          CorrStr.Xs1D = 1 - squareform(1 - CorrStr.Xs, 'tovector')';
    CorrStr.Xps = corr(X2m, 'type', 'pearson');      CorrStr.Xps1D = 1 - squareform(1 - CorrStr.Xps, 'tovector')';
    CorrStr.Xts = corr(XTm, 'type', 'pearson');      CorrStr.Xts1D = 1 - squareform(1 - CorrStr.Xts, 'tovector')';
end

function map = histmatch(pdf, xval, pdf_ref, xval_ref)
    % Parameters...
    N = numel(pdf);         L = numel(pdf_ref);
    % L = N_ref;                 if (nargin > 4),    L = varargin{1};    end
    
    % Calculate distribution mapping...
    cdf = cumDist(xval, pdf, 1);
    cdf_ref = cumDist(xval_ref, pdf_ref, 1);
    err = abs(cdf * ones(1, L) - ones(N, 1) * cdf_ref');            % N x L error matrix for matching
    [~, map] = min(err, [ ], 2);                                                  % Find best matches in CDFs
    map = (map(:) - 1) / (L - 1);                                                                % Normalization
end

function F = cumDist(x, f, varargin)
    opt = 1;    if (nargin > 2),    opt = varargin{1};      end
    switch opt
        case 1
            F = cumtrapz(x, f, 1);
            % F = diff(x(1:2)) * cumsum(f, 1);
        case 2
            % Implement simpson's rule?
        otherwise
            error("Choose other options!")
    end
end
