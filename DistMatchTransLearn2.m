function [Ypred_PS, Ypred_SS, Map, distPts, CorrStr] = DistMatchTransLearn(Xp, Yp, Xs, Ys)
    % Initial checks....
    if (size(Xp, 2) ~= size(Xs, 2)),        error('Number of features must be same!');                  end
    if (size(Xs, 1) ~= length(Ys)),         error('Secondary set must have matched datasets!');     end
    
    
    % Assign datasets...
    norm01 = @(x) (x - (ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));       % Normalize in [0, 1]
    X1 = norm01(Xp);        Y1 = (Yp(:));        n1 = size(X1, 1);
    X2 = norm01(Xs);        Y2 = (Ys(:));        n2 = size(X2, 1);
    m = size(X1, 2);
    
    tic
    %%% Transfer Learning w/ Histogram Matching...
    % Initialize parameters...
	nPts = 1e3;         distPts = linspace(0, 1, nPts+1)';                       % Sample points
    % nPts = 1e3;         distPts = struct('X', linspace(min([X1(:); X2(:)]), max([X1(:); X2(:)]) , nPts),...
    %                                                         'Y', linspace(min([Y1(:); Y2(:)]), max([Y1(:); Y2(:)]), nPts));
    L = 256;             data2idx = @(x) round((L - 1) * x) + 1;         % Convert to histogram map indices
    Hist = struct('X1', zeros(nPts, m), 'X2', zeros(nPts, m), 'X2m', zeros(nPts, m), 'X1m', zeros(nPts, m),... 
                            'Y1', zeros(nPts, 1), 'Y2', zeros(nPts, 1), 'Y2m', zeros(nPts, 1), 'Y1m', zeros(nPts, 1));
    Map = struct('X', zeros(L, m), 'Y', zeros(L, 1), 'Yinv', zeros(L, 1));
        
    %%% Covariate Map [DS1 => DS2]: per feature matching...
    X2m = zeros(n1, m);
    for j = 1 : m
        Hist.X2(:, j) = histcounts(X2(:, j), distPts, 'Normalization', 'probability');
        Hist.X1(:, j) = histcounts(X1(:, j), distPts, 'Normalization', 'probability');
        % Hist.X2(:, j) = ksdensity(X2(:, j), distPts, 'kernel', 'Normal');
        % Hist.X1(:, j) = ksdensity(X1(:, j), distPts, 'kernel', 'Normal');
        [~, Map.X(:, j)] = histeq(X1(:, j), Hist.X2(:, j));              % Mapping for j-th feature
        X2m(:, j) = Map.X(data2idx(X1(:, j)), j);                        % Mapped primary X
        X2m(X2m(:, j) < 0, j) = 0;      X2m(X2m(:, j) > 1, j) = 1;
        % Hist.X2m(:, j) = ksdensity(X2m(:, j), distPts, 'kernel', 'Normal');
        Hist.X2m(:, j) = histcounts(X2m(:, j), distPts, 'Normalization', 'probability');
    end
    
    %%% Response Maps...
    % Hist.Y2 = ksdensity(Y2, distPts, 'kernel', 'Normal');
    % Hist.Y1 = ksdensity(Y1, distPts, 'kernel', 'Normal');
    Hist.Y2 = histcounts(Y2, distPts, 'Normalization', 'probability');
    Hist.Y1 = histcounts(Y1, distPts, 'Normalization', 'probability');
    [~, Map.Y] = histeq(Y1, Hist.Y2);            Map.Y = Map.Y(:);              % Forward map: [DS1 => DS2]
    [~, Map.Yinv] = histeq(Y2, Hist.Y1);        Map.Yinv = Map.Yinv(:);     % Inverse map: [DS2 => DS1]
    Y2m = Map.Y(data2idx(Y1));                                            % Mapped primary Y
    Y1m = Map.Yinv(data2idx(Y2));                                        % Inverse mapped secondary Y
    % Hist.Y2m = ksdensity(Y2m, distPts, 'kernel', 'Normal');
    % Hist.Y1m = ksdensity(Y1m, distPts, 'kernel', 'Normal');
    Hist.Y2m = histcounts(Y2m, distPts, 'Normalization', 'probability');
    Hist.Y1m = histcounts(Y1m, distPts, 'Normalization', 'probability');
    
    %%% TL prediction...
    TrainX = (X2);	TrainY = Y2;        TestX = (X2m);
    rng(0);             nTree = 200;        RF = TreeBagger(nTree, TrainX, TrainY, 'method', 'regression');    
    PredY = predict(RF, TestX);         Y2p = PredY;        Y2p(Y2p < 0) = 0;       Y2p(Y2p > 1) = 1;
    Y1p = Map.Yinv(data2idx(Y2p));                                        % Map back to primary space
    %%% TL model finished.
    toc
    
    %%% Outputs...
    Ypred_PS = Y1p;      Ypred_SS = Y2p;
    
    clearvars CorrStr
    CorrStr.Xp = corr(X1, 'type', 'pearson');          CorrStr.Xp1D = 1 - squareform(1 - CorrStr.Xp, 'tovector')';
    CorrStr.Xs = corr(X2, 'type', 'pearson');          CorrStr.Xs1D = 1 - squareform(1 - CorrStr.Xs, 'tovector')';
    CorrStr.Xps = corr(X2m, 'type', 'pearson');      CorrStr.Xps1D = 1 - squareform(1 - CorrStr.Xps, 'tovector')';
end