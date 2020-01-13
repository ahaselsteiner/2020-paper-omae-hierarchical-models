coeffsA = [3.62, 5.77];
coeffsB = [3.54, 5.31];
coeffsC = [2.71, 6.51];

clear pi
G = 9.81;
sA = @(hs) 2 * pi * hs ./ (G * (coeffsA(1) + coeffsA(2) * sqrt(hs ./ G)).^2);
sB = @(hs) 2 * pi * hs ./ (G * (coeffsB(1) + coeffsB(2) * sqrt(hs ./ G)).^2);
sC = @(hs) 2 * pi * hs ./ (G * (coeffsC(1) + coeffsC(2) * sqrt(hs ./ G)).^2);

hs = [0:0.1:15];
fitLineWidth = 2;
markerColor = [0.3 0.3 0.3];

figSteepness = figure('position', [100 100 450 400], 'renderer', 'Painters');
hold on
load datasets-provided-ABCDEF.mat
scatter(A.Hs, computeSteepness(A.Hs, A.Tz), markerSize, markerColor, 'o')
scatter(B.Hs, computeSteepness(B.Hs, B.Tz), markerSize, markerColor, 'd')
scatter(C.Hs, computeSteepness(C.Hs, C.Tz), markerSize, markerColor, 's')
plot(hs, sA(hs), '-k', 'linewidth', fitLineWidth);
plot(hs, sB(hs), '-r', 'linewidth', fitLineWidth);
plot(hs, sC(hs), '-b', 'linewidth', fitLineWidth);
plot([0 15], [1/15 1/15], '--k');
text(13.5, 1/15, '1/15', 'verticalalignment', 'bottom');
plot([0 15], [1/30 1/30], '--k');
text(13.5, 1/30, '1/30', 'verticalalignment', 'bottom');
ylabel('steepness (-)');
xlabel('significant wave height (m)');


legend({'Dataset A', 'Dataset B', 'Dataset C', 'Model for A', 'Model for B', 'Model for C'}, ...
'location', 'southeast', 'NumColumns',2);
legend box off

function s = computeSteepness(hs, tz)
    G = 9.81;
    s = 2 * pi * hs ./ (G * tz.^2);
end