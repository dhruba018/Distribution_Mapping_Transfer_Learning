% clc;    clear;      close all

% % x = rgb2gray(imread('C:\Users\sdhruba\Downloads\Unequalized_Hawkes_Bay_NZ.jpg'));
% x = imread('pout.tif');
% % x = uint8(randi([50, 200], [8, 8]));
% % x = x/max(max(x));
% N = size(x, 1)*size(x, 2);
% minx = min(min(x));
% maxx = max(max(x));
% 
% L = 256;    Lmin = 0;   Lmax = L - 1;
% fx = zeros(L, 1);
% for k = minx:maxx
%     ind = find(x(:) == k);
%     fx(k) = length(ind)/N;
% end
% Fx = cumsum(fx);
% 
% y = uint8(zeros(size(x)));
% for i = 1:size(x, 1)
%     for j = 1:size(x, 2)
%         y(i, j) = (Lmax - Lmin)*Fx(x(i, j));
%     end
% end
% 
% figure
% subplot(121),   imshow(x),  title('Original')
% subplot(122),   imshow(y),  title('After HE')
% 
% figure
% subplot(121),   histogram(x, L),    title('Original')
% subplot(122),   histogram(y, L),    title('After HE')

% x = double(imread('pout.tif'));
% x = x/max(x(:));
% L = 256;
% % hx = histogram(x, L);
% % fx = hx.Values';
% fx = imhist(x, L)/(size(x, 1)*size(x, 2));
% Fx = cumsum(fx);
% 
% % y = histeq(x);
% % figure,     histogram(y, L)
% % y = rand(round(L/2), round(L/3));
% % y = rand(size(x));
% y = double(imread('tire.tif'));
% % y = y/max(x(:));
% y = y/max(y(:));
% yq = histeq(y, fx);
% % figure
% % hy = histogram(yq, L);
% % fy = hy.Values';
% fy = imhist(yq, L)/(size(yq, 1)*size(yq, 2));
% Fy = cumsum(fy);
% 
% figure
% subplot(121),   imshow(x)
% subplot(122),   imshow(yq)
% 
% figure
% subplot(211),   histogram(x, L)
% subplot(212),   histogram(yq, L)

% x = g1;
x = double(imread('pout.tif'));
% x = randn(1e3, 1);      x = (x - min(x))/max(x - min(x));
% x = x/max(x(:));
L = 500;
fx = imhist(x, L)/numel(x);
Fx = cumsum(fx);

% y = g2;
y = double(imread('tire.tif'));
% y = rand(200, 1);      y = (y - min(y))/max(y - min(y));
% y = y/max(x(:));
% y = y/max(y(:));
fy = imhist(y, L)/numel(y);
Fy = cumsum(fy);

% Fx = ecdf(x);
% Fy = ecdf(y);

M = zeros(size(Fy));
yq1 = zeros(size(y), 'uint8');
ind3 = cell(L, 1);
for k = 1:length(Fy)
    [~, ind] = min(abs(Fy(k) - Fx));
    M(k) = ind - 1;
%     ind3{k} = find(y == (k-1));
%     yq1(ind3{k}) = uint8(x(M(k) + 1));
end
% yq1 = M(floor(y*(L-1)) + 1) / L;
% yq1 = uint8(M(floor(y) + 1) / L);
% yq1 = uint8(x(M + 1));
yq1 = uint8(zeros(size(y)));
yq1(1:L) = uint8(x(M+1));
% yq1(ind3) = uint8(x(M+1));
yq2 = histeq(uint8(y));

% figure
% subplot(121),   imshow(x)
% subplot(122),   imshow(yq)

figure
subplot(231),       imshow(x, [ ]),          title('Original')
subplot(234),       histogram(x, L),        title('Original')
subplot(232),       imshow(yq1, [ ]),       title('Mine')
subplot(235),       histogram(yq1, L),     title('Mine')
subplot(233),       imshow(yq2, [ ]),       title('MATLAB')
subplot(236),       histogram(yq2, L),     title('MATLAB')

% figure,     histogram(x, L)
% hold on,    histogram(yq1, L)
% histogram(yq2, L),   hold off
% 
% figure
% subplot(211),	histogram(x, L)
% subplot(212),   histogram(yq1, L),  hold on
% histogram(yq2, L),   hold off


%%
clear
% y = im2double(imread('tire.tif'));
y = im2double(imread('pout.tif'));
[m, n] = size(y);
L = 100;

figure(1)
subplot(121),	imshow(y, [ ])
subplot(122),	histogram(y, L, 'Normalization', 'probability')

t = rand(m, n);
figure(2)
subplot(121),	imshow(t, [ ])
subplot(122),	histogram(t, L, 'Normalization', 'probability')

fy = imhist(y, L) / (m*n);
Fy = cumsum(fy);

fy2 = histcounts(y(:), L)' / (m*n);
Fy2 = cumsum(fy2);

figure(3)
subplot(121),   bar(linspace(0, 1, L), fy)
subplot(122),   bar(linspace(0, 1, L), fy2)

ft2 = histcounts(t(:), L)' / (m*n);
Ft2 = cumsum(ft2);

figure(4)
subplot(121),   bar(linspace(0, 1, L), fy2)
subplot(122),   bar(linspace(0, 1, L), ft2)

M = zeros(L, 1);
lvls = linspace(0, 1, L)';
Th = 1e-2;
yq = zeros(m, n);
for k = 1 : L
    [~, ind] = min(abs(fy2(k) - ft2));
    M(k) = ind - 1;
    idx1 = find(abs(y - lvls(k)) <= Th);
%     idx2 = find(abs(y - lvls(k)) <= Th);
%     yq(idx1) = lvls(ind);
end
% yq = interp
% yq = yq(:);
% Idx = find(yq ~= 0);
% Idxq = 1:(m*n);
% yq = interp1(Idx, yq(yq ~= 0), Idxq);
% yq = knnimpute(reshape(yq, [m, n]), 10);
% Idx1 = [ ];     [Idx1(:, 1), Idx1(:, 2)] = find(yq ~= 0);
% yq = interp2(Idx1, yq(Idx1), [1:m, 1:n]);
yq = M(floor(y) + 1);
fq2 = histcounts(yq(:), L)' / (m*n);

yh = histeq(y, ft2);
fh2 = histcounts(yh(:), L)' / (m*n);

figure(5)
subplot(141),   bar(lvls, fy2)
subplot(142),   bar(lvls, ft2)
subplot(143),   bar(lvls, fq2)
subplot(144),   bar(lvls, fh2)

figure(6)
subplot(141),   imshow(y, [ ])
subplot(142),   imshow(t, [ ])
subplot(143),   imshow(yq, [ ])
subplot(144),   imshow(yh, [ ])

%%
x1 = 0.2 + 0.03*randn(1e3, 1);
% x2 = -0.2 + 0.4*rand(1e3, 1);

x2 = 0.2 + 0.03*rand(1e3, 1);
% - 0.4 + 0.05*rand(1e3, 1);

N = 250;
x1q = histeq(x1, histcounts(x2, N));

figure(7)
subplot(131),   histogram(x1, N)
subplot(132),   histogram(x2, N)
subplot(133),   histogram(x1q, N)

f1 = histcounts(x1, N) / numel(x1);
f2 = histcounts(x2, N) / numel(x2);
F1 = cumsum(f1);
F2 = cumsum(f2);
M = zeros(size(F1));
for k = 1:N
    [~, ind] = min(abs(F1(k) - F2));
    M(k) = ind - 1;
end
x1q2 = M(floor(x1*(N-1)) + 1) / N;
figure(7)
subplot(141),   histogram(x1, N)
subplot(142),   histogram(x2, N)
subplot(143),   histogram(x1q, N)
subplot(144),   histogram(x1q2, N)

%%
nBin = 50;
figure(1),
subplot(121),       histogram(dr3, nBin,  'Normalization', 'probability'),      title('dr3')
subplot(122),       histogram(dr2, nBin,  'Normalization', 'probability'),      title('dr2')

[fd3, v3] = histcounts(dr3, nBin);      fd3 = fd3' / numel(dr3);    v3 = cumsum(diff(v3))';
[fd2, v2] = histcounts(dr2, nBin);      fd2 = fd2' / numel(dr2);    v2 = cumsum(diff(v2))';

figure(2),
subplot(121),       bar(v3, fd3),      title('dr3 - bar')
subplot(122),       bar(v2, fd2),      title('dr2 - bar')

[dr3q, T] = histeq(dr3, fd2);
[fd3q, v3q] = histcounts(dr3, nBin);      fd3q = fd3q' / numel(dr3);    v3q = cumsum(diff(v3q))';

figure(3),
subplot(131),       histogram(dr3, nBin,  'Normalization', 'probability'),        title('dr3')
subplot(132),       histogram(dr2, nBin,  'Normalization', 'probability'),        title('dr2')
subplot(133),       histogram(dr3q, nBin,  'Normalization', 'probability'),      title('dr3q')



%%
clc;    clear;      clf

L = 256;    xBins = linspace(0, 1, L+1)';

x1 = im2double(imread('pout.tif'));
fx1 = histcounts(x1, xBins, 'Normalization', 'probability')';

x2 = im2double(imread('tire.tif'));
fx2 = histcounts(x2, xBins, 'Normalization', 'probability')';

x1m = histeq(x1, fx2);
fx1m = histcounts(x1m, xBins, 'Normalization', 'probability')';

x1m2 = HistMatch(x1, fx2, xBins);
fx1m2 = histcounts(x1m2, xBins, 'Normalization', 'probability')';

figure(1),          binColorRGB = [2, 6, 6] / 10;
subplot(2,4,1),     imshow(x1, [ ]),    title('Target')
subplot(2,4,5),     bar(xBins(1:end-1), fx1, 'FaceColor', binColorRGB),     xlim([0, 1])
subplot(2,4,2),     imshow(x2, [ ]),    title('Reference')
subplot(2,4,6),     bar(xBins(1:end-1), fx2, 'FaceColor', binColorRGB),     xlim([0, 1])
subplot(2,4,3),     imshow(x1m, [ ]),   title('Matched')
subplot(2,4,7),     bar(xBins(1:end-1), fx1m, 'FaceColor', binColorRGB),    xlim([0, 1])
subplot(2,4,4),     imshow(x1m2, [ ]),	title('Matched2')
subplot(2,4,8),     bar(xBins(1:end-1), fx1m2, 'FaceColor', binColorRGB),	xlim([0, 1])


function matched = HistMatch(target, refDist, varargin)
    % Parameters...
    if (nargin > 2),    bins = varargin{1}(:);     L = numel(bins) - 1;
    else,               L = 256;                   bins = linspace(0, 1, L+1)';     end
    
    % Histogram matching...
    targetDist = histcounts(target, bins, 'Normalization', 'probability')';
    targetCumDist = cumsum(targetDist);     refCumDist = cumsum(refDist);    
    map = zeros(L);
    for k = 1 : L
        [~, ind] = min(abs(targetCumDist(k) - refCumDist));
        map(k) = ind - 1;
    end
    matched = map(floor(target * (L - 1)) + 1) / L;
end




















