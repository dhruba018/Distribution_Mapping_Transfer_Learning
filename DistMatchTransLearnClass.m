function [Ypred, Map, distPts, CorrStr] = DistMatchTransLearnClass(Xp, Yp, Xs, Ys, mdlChoice)
    % Initial checks....
    if (size(Xp, 2) ~= size(Xs, 2)),        error('Number of features must be same!');                  end
    if (size(Xs, 1) ~= length(Ys)),         error('Secondary set must have matched datasets!');     end
    
    % Assign datasets...
    norm01 = @(x) (x - min(ones(size(x, 1), 1) * min(x, [ ], 1))) ./ (ones(size(x, 1), 1) * range(x, 1));       % Normalize in [0, 1]
    X1 = norm01(Xp);        Y1 = (Yp(:));        n1 = size(X1, 1);
    X2 = norm01(Xs);        Y2 = (Ys(:));        n2 = size(X2, 1);
    m = size(X1, 2);
    
    tic
    %%% Transfer Learning w/ Histogram Matching...
    % Initialize parameters...
	nPts = 0.5e3;         distPts = linspace(0, 1, nPts)';                       % Sample points
    L = 256;                data2idx = @(x) round((L - 1) * x) + 1;          % Convert to histogram map indices
    Hist = struct('X1', zeros(nPts, m), 'X2', zeros(nPts, m), 'X2m', zeros(nPts, m), 'X1m', zeros(nPts, m));
    Map = struct('X', zeros(L, m));
        
    %%% Covariate Map [DS1 => DS2]: per feature matching...
    X2m = zeros(n1, m);
    for j = 1 : m
        Hist.X2(:, j) = ksdensity(X2(:, j), distPts, 'kernel', 'Normal');
        Hist.X1(:, j) = ksdensity(X1(:, j), distPts, 'kernel', 'Normal');
        [~, Map.X(:, j)] = histeq(X1(:, j), Hist.X2(:, j));              % Mapping for j-th feature
        X2m(:, j) = Map.X(data2idx(X1(:, j)), j);                        % Mapped primary X
        Hist.X2m(:, j) = ksdensity(X2m(:, j), distPts, 'kernel', 'Normal');
    end
    
    %%% TL prediction...
    TrainX = zscore(X2);            TrainY = Y2;            TestX = zscore(X2m);
%     mdlChoice = 3;
    switch mdlChoice
        case 1
            rng(0);                     nTree = 100;            RF = TreeBagger(nTree, TrainX, TrainY, 'method', 'classification');
            PredY = predict(RF, TestX);                     Y1p = str2num(cell2mat(PredY));
        case 2
            KNN = fitcknn(TrainX, TrainY, 'NumNeighbors', 7);
            PredY = predict(KNN, TestX);                  Y1p = PredY;
        case 3
            Kernel = 'poly';        PolyOrd = 3;
            SVM = fitcsvm(TrainX, TrainY, 'KernelFunction', Kernel, 'PolynomialOrder', PolyOrd);
            PredY = predict(SVM, TestX);                  Y1p = PredY;
    end
    
    %%% TL model finished.
    toc
    
    %%% Outputs...
    Ypred = Y1p;
    
    clearvars CorrStr
    CorrStr.Xp = corr(X1, 'type', 'pearson');          CorrStr.Xp1D = 1 - squareform(1 - CorrStr.Xp, 'tovector')';
    CorrStr.Xs = corr(X2, 'type', 'pearson');          CorrStr.Xs1D = 1 - squareform(1 - CorrStr.Xs, 'tovector')';
    CorrStr.Xps = corr(X2m, 'type', 'pearson');      CorrStr.Xps1D = 1 - squareform(1 - CorrStr.Xps, 'tovector')';
    