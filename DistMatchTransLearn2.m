function [Ypred_PS, Ypred_SS, Map, distPts, CorrStr] = DistMatchTransLearn2(Xp, Yp, Xs, Ys)
    % Initial checks....
    if (size(Xp, 2) ~= size(Xs, 2)),        error('Number of features must be same!');                  end
    if (size(Xs, 1) ~= length(Ys)),         error('Secondary set must have matched datasets!');     end
    
    % Assign datasets...
    norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));       % Normalize in [0, 1]
    X1 = norm01(Xp);        Y1 = (Yp(:));        n1 = size(X1, 1);
    X2 = norm01(Xs);        Y2 = (Ys(:));        n2 = size(X2, 1);
    m = size(X1, 2);
    
    % nPts = 1e3;         Y_distPts = [linspace(min(Y1), max(Y1), nPts); linspace(min(Y2), max(Y2), nPts)]';
    
    tic
    %%% Transfer Learning w/ Histogram Matching...
    % Initialize parameters...
	nPts = 1e3;         distPts = linspace(0, 1, nPts)';                       % Sample points
    % nPts = 1e3;         distPts = struct('X', linspace(min([X1(:); X2(:)]), max([X1(:); X2(:)]) , nPts),...
    %                                                         'Y', linspace(min([Y1(:); Y2(:)]), max([Y1(:); Y2(:)]), nPts));
    MapCalc = 1;      L = 256 * (MapCalc == 0) + nPts * (MapCalc ~= 0);
    data2idx = @(z) round((L - 1) * z) + 1;         % Convert to histogram map indices
    Hist = struct('X1', zeros(nPts, m), 'X2', zeros(nPts, m), 'X2m', zeros(nPts, m), 'X1m', zeros(nPts, m),... 
                            'Y1', zeros(nPts, 1), 'Y2', zeros(nPts, 1), 'Y2m', zeros(nPts, 1), 'Y1m', zeros(nPts, 1));
    Map = struct('X', zeros(L, m), 'Y', zeros(L, 1), 'Yinv', zeros(L, 1));
    
    %%% Covariate Map [DS1 => DS2]: per feature matching...
    X2m = zeros(n1, m);
    for j = 1 : m
        Hist.X2(:, j) = ksdensity(X2(:, j), distPts, 'kernel', 'Normal');
        Hist.X1(:, j) = ksdensity(X1(:, j), distPts, 'kernel', 'Normal');
        
        % Mapping for j-th feature...
        switch MapCalc
            case 0
                [~, Map.X(:, j)] = histeq(X1(:, j), Hist.X2(:, j));
            case 1
                Map.X(:, j) = histmatch(Hist.X1(:, j), distPts, Hist.X2(:, j), distPts);
        end
        X2m(:, j) = Map.X(data2idx(X1(:, j)), j);                        % Mapped primary X
        Hist.X2m(:, j) = ksdensity(X2m(:, j), distPts, 'kernel', 'Normal');
    end
    
    %%% Response Maps...
    Hist.Y2 = ksdensity(Y2, distPts, 'kernel', 'Normal');
    Hist.Y1 = ksdensity(Y1, distPts, 'kernel', 'Normal');
    
    % Forward & Inverse Maps...
    switch MapCalc
        case 0
            [~, Map.Y] = histeq(Y1, Hist.Y2);            Map.Y = Map.Y(:);              % DS1 => DS2
            [~, Map.Yinv] = histeq(Y2, Hist.Y1);        Map.Yinv = Map.Yinv(:);     % DS2 => DS1
        case 1
            Map.Y = histmatch(Hist.Y1, distPts, Hist.Y2, distPts);                      % DS1 => DS2
            Map.Yinv = histmatch(Hist.Y2, distPts, Hist.Y1, distPts);                  % DS2 => DS1
    end
    Y2m = Map.Y(data2idx(Y1));                                            % Mapped primary Y
    Y1m = Map.Yinv(data2idx(Y2));                                        % Inverse mapped secondary Y
    Hist.Y2m = ksdensity(Y2m, distPts, 'kernel', 'Normal');
    Hist.Y1m = ksdensity(Y1m, distPts, 'kernel', 'Normal');
    
    %%% TL prediction...
    TrainX = X2;	TrainY = Y2;        TestX = X2m;
    rng(0);             nTree = 200;        RF = TreeBagger(nTree, TrainX, TrainY, 'method', 'regression');    
    PredY = predict(RF, TestX);         Y2p = PredY;       % Y2p(Y2p < 0) = 0;       Y2p(Y2p > 1) = 1;
    Y1p = Map.Yinv(data2idx(Y2p));                                        % Map back to primary space
    % Y1p = Map.
    %%% TL model finished.
    toc
    
    %%% Outputs...
    Ypred_PS = Y1p;      Ypred_SS = Y2p;
    
    clearvars CorrStr
    CorrStr.Xp = corr(X1, 'type', 'pearson');          CorrStr.Xp1D = 1 - squareform(1 - CorrStr.Xp, 'tovector')';
    CorrStr.Xs = corr(X2, 'type', 'pearson');          CorrStr.Xs1D = 1 - squareform(1 - CorrStr.Xs, 'tovector')';
    CorrStr.Xps = corr(X2m, 'type', 'pearson');      CorrStr.Xps1D = 1 - squareform(1 - CorrStr.Xps, 'tovector')';
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
    map = (map(:) - 1) / (L - 1);                                                 % Normalization
end

function F = cumDist(x, f, varargin)
    opt = 1;    if (nargin > 2),    opt = varargin{1};      end
    switch opt
        case 1
            F = cumtrapz(x, f, 1);
        case 2
            % Implement simpson's rule?
        otherwise
            error("Choose other options!")
    end
end
