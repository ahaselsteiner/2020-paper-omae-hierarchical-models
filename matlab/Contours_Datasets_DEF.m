addpath('C:\Users\Andreas Haselsteiner\Documents\MATLAB\2019-paper-predicting-wave-heights\work-spaces');
addpath('C:\Users\Andreas Haselsteiner\Documents\MATLAB\2019-paper-predicting-wave-heights\exponentiated-weibull');
addpath('C:\Users\Andreas Haselsteiner\Documents\MATLAB\compute-hdc');
load datasets-provided-ABCDEF.mat
DELTA = 5;
DATASET_CHAR = 'F';


binWidth = 2;
minDataPointsInBin = 20;

if DATASET_CHAR == 'D'
    v = D.V;
    hs = D.Hs;
elseif DATASET_CHAR == 'E'
    v = E.V;
    hs = E.Hs;
elseif DATASET_CHAR == 'F'
    v = F.V;
    hs = F.Hs;
end


%Is it Weibull distributed or do we need a exponent? For dataset F: Yes
pdV_ExpWbl = ExponentiatedWeibull();
pdV_ExpWbl.fitDist(v, 'WLS')

xi = sort(v);
n = length(v);
i = [1:1:n]';
pi = (i - 0.5) / n;
pstar_i = computePStar(pi, 1);
xHat_ExpMLE = pdV_ExpWbl.icdf(pi);

figExpWbl = figure('position', [100 100 450 400]);
% 2 parameter Weibull
parmHat = wblfit(v);
pd2PWbl = ExponentiatedWeibull(parmHat(1), parmHat(2), 1);
xHat_2pWbl = pd2PWbl.icdf(pi);
plot(log10(xi), pstar_i, 'xk', 'markersize', 4, 'color', [0.4 0.4 0.4]);
hold on
plot(log10(xHat_2pWbl), pstar_i, '--k', 'linewidth', 1.5);
plot(log10(xHat_ExpMLE), pstar_i, '-b', 'linewidth', 1.5);
xticks = [0.5 1 2 4 8 12 20 30 40];
set(gca, 'xtick', log10(xticks));
set(gca, 'xticklabel', xticks);
xlabel('wind speed (m/s)');
yticks =   [0.01 0.05 0.1 0.2 0.5 0.8 0.9 0.95 0.99 0.999 0.9999 0.99999 0.999999];
set(gca, 'ytick', computePStar(yticks, 1));
set(gca, 'yticklabel', yticks);
ylabel('probability (-)');
xlim(log10([2 35]));
ylim(computePStar([0.1 0.9999999], 1))
grid on
format3Digits = '%.3f';
twoParamString = ['2-parameter Weibull ' char(10) '(\alpha = '...
    f3digits(parmHat(1)) ', \beta = ' f3digits(parmHat(2)) ')'];
expString = ['Exponentiated Weibull ' char(10) ...
    '(\alpha = ' f3digits(pdV_ExpWbl.Alpha) ...
    ', \beta = ' f3digits(pdV_ExpWbl.Beta) ...
    ', \delta = ' f3digits(pdV_ExpWbl.Delta) ')'];
legend({['Dataset ' DATASET_CHAR], twoParamString, expString}, 'location', 'northwest');
legend box off



nOfBins = ceil((max(v) / binWidth));

pdHs_Exp = ExponentiatedWeibull.empty(nOfBins, 0);
xi = cell(nOfBins, 1);
pstar_i = cell(nOfBins, 1);
xHat_2pWbl = cell(nOfBins, 1);
xHat_Exp = cell(nOfBins, 1);
hsInBin = cell(nOfBins,1);
lowerLimit = nan(nOfBins, 1);
upperLimit = nan(nOfBins, 1);
medianOfHs = nan(nOfBins, 1);
varOfHs = nan(nOfBins, 1);



for i = 1:nOfBins
    lowerLimit(i) = (i - 1) * binWidth;
    upperLimit(i) = i * binWidth;
    vIsInBin = (v > lowerLimit(i)) .* (v <= upperLimit(i));
    hsInBin{i} = hs(logical(vIsInBin));

    if length(hsInBin{i}) >= minDataPointsInBin
        xi{i} = sort(hsInBin{i});
        n = length(xi{i});
        j = [1:1:n]';
        pi = (j - 0.5) / n;
        pstar_i{i} = computePStar(pi, 1);
        
        pdHs_Exp(i) = ExponentiatedWeibull();
        pdHs_Exp(i).fitDist(hsInBin{i}, 'WLS', 'delta', DELTA);
        parmHat = wblfit(hsInBin{i});
        pdHs_2pWbl(i) = makedist('Weibull', parmHat(1), parmHat(2));

        xHat_2pWbl{i} = pdHs_2pWbl(i).icdf(pi);
        xHat_Exp{i} = pdHs_Exp(i).icdf(pi);
        
        % Compute some metrics for analysis
        medianOfHs(i) = median(hsInBin{i});
        varOfHs(i) = var(hsInBin{i});
    end
end

medianOfHs = medianOfHs(~isnan(medianOfHs));
varOfHs = varOfHs(~isnan(varOfHs));

alphas = nan(length(pdHs_Exp), 1);
betas = nan(length(pdHs_Exp), 1);
deltas = nan(length(pdHs_Exp), 1);
binCenter = nan(length(pdHs_Exp), 1);
for i = 1:length(pdHs_Exp)
    alphas(i) = pdHs_Exp(i).Alpha;
    betas(i) = pdHs_Exp(i).Beta;
    deltas(i) = pdHs_Exp(i).Delta;
    binCenter(i) = (i - 0.5) * binWidth;
end

figMarginalHsFits = figure('position', [100, 100, 1400, 400], 'renderer', 'Painters');
doPlot = 0;
for i = 1:length(pdHs_Exp)
        if binCenter(i) == 1 || binCenter(i) == 11 || binCenter(i) == 21
            doPlot = 1;
        else
            doPlot = 0;
        end
        if binCenter(i) == 1
            plotI = 1;
            my_xlim = [0.1 10];
        elseif binCenter(i) == 11
            plotI = 2;
            my_xlim = [1 11];
        elseif binCenter(i) == 21
            plotI = 3;
            my_xlim = [3 15];
        end
        if doPlot
            subplot(1, 3, plotI);
            xlim(log10(my_xlim));
            plot(log10(xi{i}), pstar_i{i}, 'xk', 'markersize', 4, 'color', [0.3 0.3 0.3]);
            hold on
            plot(log10(xHat_2pWbl{i}), pstar_i{i}, '--k', 'linewidth', 1.5);
            plot(log10(xHat_Exp{i}), pstar_i{i}, '-b', 'linewidth', 1.5);
            xticks = [0.1 0.2 0.5 1 2 4 6 8 10 12];
            set(gca, 'xtick', log10(xticks));
            set(gca, 'xticklabel', xticks);
            yticks = [0.05 0.1 0.2 0.5 0.9 0.99 0.9999];
            set(gca, 'ytick', computePStar(yticks, 1));
            set(gca, 'yticklabel', yticks);
            xlabel('significant wave height (m)');
            if plotI == 1 
                ylabel('probability (-)');
            end
            ylim(computePStar([0.05 0.999999], 1))
            grid on
            twoParamString = ['2-parameter Weibull ' char(10) '(\alpha = '...
                f3digits(pdHs_2pWbl(i).A) ', \beta = ' f3digits(pdHs_2pWbl(i).B) ')'];
            expString = ['Exponentiated Weibull ' char(10) ...
                '(\alpha = ' f3digits(pdHs_Exp(i).Alpha) ...
                ', \beta = ' f3digits(pdHs_Exp(i).Beta) ...
                ', \delta = ' f3digits(pdHs_Exp(i).Delta) ')'];
            legend({['Dataset ' DATASET_CHAR], twoParamString, expString}, 'location', 'southeast');
            legend box off
            title([num2str(lowerLimit(i)) ' $< v \le$' num2str(upperLimit(i)) ' m/s'], 'interpreter', 'latex');
        end
end

powerFunction = 'a + b*x^c';
logisticFunction4p = 'y0 + L / (1 + exp(-k * (x - x0)))'; % Paramterization as in https://en.wikipedia.org/wiki/Logistic_function

fo3p = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 0, 0],...
               'Upper',[Inf, Inf, Inf],...
               'StartPoint',[1 1 1]);
fo4p = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 0, 0, 0],...
               'Upper',[Inf, Inf, Inf Inf],...
               'StartPoint', [1 1 1 1]);

fBeta = fit(binCenter, betas, logisticFunction4p, fo4p);
fMedian = fit(binCenter, medianOfHs, powerFunction, fo3p);
medianCoeff = -1 * log(1 - 0.5^(1 / DELTA));
fAlpha = @(x1)fMedian(x1) ./ medianCoeff.^(1 ./ fBeta(x1));

x = [0:0.1:max(binCenter)]';
fitting_line_width = 2;
marker_face_color = [0.8 0.8 0.8];

fig = figure('position', [100 100 1100 300]);
subplot(1,3,1);
plot(binCenter, medianOfHs, 'ok', 'markerfacecolor', marker_face_color);
hold on
plot(x, fMedian(x), '-b', 'linewidth', fitting_line_width)
eq_string = ['' f3digits(fMedian.a)  ' + ' f3digits(fMedian.b) ' * x^{'  f3digits(fMedian.c) '}']
legend({'from marginal distribution', eq_string})
legend box off
xlabel('wind speed (m/s)');
ylabel('h_{s50} (m)');
box off
xlim([0 30])
ylim([0 14])



subplot(1,3,2);
plot(binCenter, betas, 'ok');
hold on
plot(x, fBeta(x), '-b', 'linewidth', fitting_line_width);
eq_string = ['' f3digits(fBeta.y0)  ' + ' f3digits(fBeta.L) ' / [1 + e^{-'  f3digits(fBeta.k) ' * (x - ' f3digits(fBeta.x0) ')}]'];
legend({'from marginal distribution', eq_string}, 'fontsize', 8)
legend box off
xlabel('wind speed (m/s)');
ylabel('\beta_{hs}');
box off
xlim([0 30])
ylim([0 3.5])

subplot(1,3,3);
plot(binCenter, alphas, 'ok');
hold on
plot(x, fAlpha(x), '-b', 'linewidth', fitting_line_width);
eq_string = 'h_{s50} / 2.0445^{(1 / \beta_{hs})}'
legend({'from marginal distribution', eq_string})
legend box off
xlabel('wind speed (m/s)');
ylabel('\alpha_{hs}');
box off
xlim([0 30])
ylim([0 10])
sgtitle(['Dataset ' DATASET_CHAR])




PM.name = 'V-Hs model';
PM.modelType = 'CMA';
PM.distributions = {'exponentiated-weibull'; 'exponentiated-weibull'};
PM.isConditionals = {[0 0 0 ]; [1 1 1]};
PM.coeffs = {
    {pdV_ExpWbl.Alpha pdV_ExpWbl.Beta pdV_ExpWbl.Delta}; 
    { 
    fAlpha;
    @(x1)fBeta(x1);
    @(x1)deltas(1)}
    };
PM(1).labels = {'Wind speed (m/s)';
    'sign. wave height (m)'};
PM.gridCenterPoints = {0.05:0.1:50; 0.05:0.1:30};

n_years = 1;
alpha = 1 / (n_years * 365.25 * 24);
[fm, x1_1Yr, x2_1Yr] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);

figContour = figure('position', [100 100 350 300]);
plot(v, hs, '.k');
hold on
plot(x1_1Yr{1}, x2_1Yr{1}, '--b', 'linewidth', 2);


n_years = 50;
alpha = 1 / (n_years * 365.25 * 24);
[fm, x1_50Yr, x2_50Yr] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);
figure(figContour);
plot(x1_50Yr{1}, x2_50Yr{1}, '-b', 'linewidth', 2);

xlabel('wind speed (m/s)');
ylabel('significant wave height (m)');

legend({'Observations', '1-year contour', '50-year contour', '500-year contour'}, ...
    'location', 'northwest');
legend box off
box off
title(['Dataset ' DATASET_CHAR])
xlim([0 35]);
ylim([0 20]);

function pstar_i = computePStar(pi, delta)
    pstar_i = log10(-log(1 - pi.^(1 ./ delta)));
end

function s = f3digits(x)
 s = sprintf('%.3g', round(x, 3, 'significant'));
end
