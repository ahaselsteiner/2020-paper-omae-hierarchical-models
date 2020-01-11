addpath('C:\Users\Andreas Haselsteiner\Documents\MATLAB\2019-paper-predicting-wave-heights\work-spaces');
addpath('C:\Users\Andreas Haselsteiner\Documents\MATLAB\2019-paper-predicting-wave-heights\exponentiated-weibull');
addpath('C:\Users\Andreas Haselsteiner\Documents\MATLAB\compute-hdc');
load datasets-provided-ABCDEF.mat

hs = C.Hs;
tz = C.Tz;



% Wave height
pdV_ExpWbl = ExponentiatedWeibull();
pdV_ExpWbl.fitDist(hs, 'WLS')
pdV_ExpWbl.qqplot(hs);

binWidth = 0.5;
minDataPointsInBin = 20;
nOfBins = ceil((max(hs) / binWidth));


pdTz = makedist('Lognormal');
pdS = makedist('Lognormal');
xi = cell(nOfBins, 1);
pstar_i = cell(nOfBins, 1);
xHat_2pWbl = cell(nOfBins, 1);
xHat_Exp = cell(nOfBins, 1);
tzInBin = cell(nOfBins,1);
steepnessInBin = cell(nOfBins,1);
lowerLimit = nan(nOfBins, 1);
upperLimit = nan(nOfBins, 1);
medianOfTz = nan(nOfBins, 1);
medianOfSteepness = nan(nOfBins, 1);
meanOfHs = nan(nOfBins, 1);
varOfTz = nan(nOfBins, 1);


fig = figure('position', [100, 100, 1400, 350], 'renderer', 'Painters');
for i = 1:nOfBins
    lowerLimit(i) = (i - 1) * binWidth;
    upperLimit(i) = i * binWidth;
    hsIsInBin = (hs > lowerLimit(i)) .* (hs < upperLimit(i));
    hsInBin{i} = hs(logical(hsIsInBin));
    tzInBin{i} = tz(logical(hsIsInBin));
    steepnessInBin{i} = 2 * pi * hsInBin{i} ./ (9.81 * tzInBin{i}.^2);

    
    if length(tzInBin{i}) >= minDataPointsInBin
        pdTz(i) = fitdist(tzInBin{i}, 'lognormal');
        pdS(i) = fitdist(steepnessInBin{i}, 'lognormal');
        
        % Compute some metrics for analysis
        medianOfTz(i) = median(tzInBin{i});
        
        medianOfSteepness(i) = median(steepnessInBin{i});
        varOfTz(i) = var(tzInBin{i});
        meanOfHs(i) = mean(hsInBin{i});
    end
end

medianOfTz = medianOfTz(~isnan(medianOfTz));
medianOfSteepness = medianOfSteepness(~isnan(medianOfSteepness));
varOfTz = varOfTz(~isnan(varOfTz));
meanOfHs = meanOfHs(~isnan(meanOfHs));

for i = 1:length(pdTz)
        subplot(2, ceil(length(pdTz) / 2), i);
        histfit(steepnessInBin{i}, 10, 'lognormal')
        title([num2str(lowerLimit(i)) ' < hs < ' num2str(upperLimit(i)) ' m']);
end


tz_mus = nan(length(pdTz), 1);
tz_sigmas = nan(length(pdTz), 1);
s_mus = nan(length(pdTz), 1);
s_sigmas = nan(length(pdTz), 1);
binCenter = nan(length(pdTz), 1);
for i = 1:length(pdTz)
    tz_mus(i) = pdTz(i).mu;
    tz_sigmas(i) = pdTz(i).sigma;
    s_mus(i) = pdS(i).mu;
    s_sigmas(i) = pdS(i).sigma;
    binCenter(i) = (i - 0.5) * binWidth;
end

fixedLinearFunction = [num2str(min(tz_sigmas)) ' + a * x'];
linearFunction = 'a + b * x';
quadraticFunction = 'a + b * x^2';
polynom2 = 'a + b * x + c*x^2';
inhsQuadraticFunction = 'a - (x + b)^(-2)';
cubicFunction = 'a + b * x^3';
powerFunction = 'a + b*x^c';
lnSquareFunction = 'log(a + b * sqrt(x / 9.81))'; % Idea based on IEC 61400-3-1:2019 p. 37
lnPowerFunction = 'log(a + b*x^c)';
%powerFunctionWMinimum = [num2str(min(tz_sigmas)) '+ b*x^c'];
powerFunctionWMinimum = ['(' num2str(min(tz_sigmas)) ' < a + b*x^c) * (a + b*x^c - ' num2str(min(tz_sigmas)) ') + ' num2str(min(tz_sigmas))];

expFunction2p = 'a * exp(b * x)';
expFunction3p = 'a + b * exp(c * x)';
logisticFunction = 'c / (1 + a * exp(-1 * r * x))'; % See https://www.classzone.com/eservices/home/pdf/student/LA208HAD.pdf
logisticFunctionWOffset = [num2str(min(tz_sigmas)) ' + a / (1 + exp(b*(-1*(x - c))))'];
logisticFunction4p = 'y0 + L / (1 + exp(-k * (x - x0)))'; % Paramterization as in https://en.wikipedia.org/wiki/Logistic_function
errorFunction = 'a + b * erf(x * c)';
steepnessFunction =   'log(sqrt(2 * pi * x / (9.81 * s)))';
steepnessFunction2p = 'log(sqrt(2 * pi * x / (9.81 * (a * log(x + b)))))';
steepnessFunction3p = 'log(sqrt(2 * pi * x / (9.81 * (a + b * x^c))))';
%steepnessFunction3p = 'log(sqrt(2 * pi * x / (9.81 * (b * (x + a)^c))))';
steepnessLogistic = 'log(sqrt(2 * pi * x / (9.81 * (c / (1 + a * exp(-1 * r * x))))))';

fo1p = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'StartPoint',[1]);
fo2p = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, -Inf],...
               'Upper',[Inf, Inf],...
               'StartPoint',[1 1], ...
               'weights', 1./ medianOfTz);
fo2pSpecial = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 1],...
               'Upper',[100, 100],...
               'StartPoint',[1 1], ...
               'weights', 1./ medianOfTz);
fo3p = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 0, 0],...
               'Upper',[Inf, Inf, Inf],...
               'StartPoint',[1 1 1], ...
               'weights', 1./ medianOfTz);
fo3pSpecial = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0, 0, -100],...
            'Upper',[100, 100, 100],...
            'StartPoint',[1 1 1], ...
            'weights', 1./ medianOfTz);
fo3pSteepness = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0, 0, -100],...
            'Upper',[100, 100, 100],...
            'StartPoint',[1 1 1], ...
            'weights', 1./ medianOfSteepness);
fo3pAllowNegative = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf, -Inf, -Inf],...
               'Upper',[Inf, Inf, Inf],...
            'StartPoint',[1 1 1], ...
            'weights', 1./ medianOfTz);
fo4p = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf, -Inf, -Inf, -Inf],...
               'Upper',[Inf, Inf, Inf Inf],...
               'StartPoint', [1 1 1 1], ...
               'weights', 1 ./ tz_sigmas);

fMuTz = fit(meanOfHs, tz_mus, lnSquareFunction, fo2p)
%fMuTz = fit(meanOfHs, tz_mus, steepnessFunction3p, fo3p);
%fMuTz = fit(meanOfHs, tz_mus, steepnessFunction2p, fo2pSpecial);
%fMuTz = fit(meanOfHs, tz_mus, powerFunction, fo3pSpecial);
%myFunc = 'a * log(x + b)';
%fMuTz = fit(meanOfHs, tz_mus, myFunc, fo2pSpecial);

%fMedian = fit(meanOfHs, medianOfSteepness, myFunc, fo2p); 
%fPositiveMedian = @(x) double(fMedian(x) > 0) .* fMedian(x)
%fMuTz = @(x) log(sqrt(2 * pi * x ./ (9.81 * fPositiveMedian(x)')));

powerdecrease3 = 'a + 1 / ((x + b)^c)'
asymdecrease3 = 'a + 1 / (b * (x + c))'
fSigmaTz = fit(meanOfHs, tz_sigmas, asymdecrease3, fo3p);


% fMuS = fit(meanOfHs, s_mus, lnPowerFunction, fo3pSteepness);
% fSigmaS = fit(meanOfHs, s_sigmas, expFunction3p, fo3p);

steepness = 2 * pi * hs ./ (9.81 * tz.^2);
figSteepness = figure()
plot(hs, steepness, '.k');;
hold on
plot(meanOfHs, medianOfSteepness, 'or');
x = [0:0.1:max(hs)];

xlabel('hs');
ylabel('steepness');


figure
subplot(2,1,1);
plot(meanOfHs, tz_mus, 'ok');
hold on
plot(fMuTz);
%x = [0: 0.1: max(hs)];
%plot(x, fMuTz(x))
ylabel('mu');

subplot(2,1,2);
plot(meanOfHs, tz_sigmas, 'ok');
hold on
plot(fSigmaTz);
ylabel('sigma');
xlabel('Signifi speed (m/s)');

% figure
% subplot(2,1,1);
% plot(meanOfHs, s_mus, 'ok');
% hold on
% plot(fMuS);
% ylabel('mu');
% 
% subplot(2,1,2);
% plot(meanOfHs, s_sigmas, 'ok');
% hold on
% plot(fSigmaS);
% ylabel('sigma');
% xlabel('Wind speed (m/s)');

PM.name = 'Hs-Tz model';
PM.modelType = 'CMA';
PM.distributions = {'exponentiated-weibull'; 'lognormal'};
PM.isConditionals = {[0 0 0 ]; [1 1]};
PM.coeffs = {
    {pdV_ExpWbl.Alpha pdV_ExpWbl.Beta pdV_ExpWbl.Delta}; 
    { 
    @(x1)fMuTz(x1);
    @(x1)fSigmaTz(x1)}
    };
PM(1).labels = {'Sign. wave height (m)';
    'Zero-upcrossing period (s)'};
PM.gridCenterPoints = {0.05:0.1:50; 0.05:0.1:30};

n_years = 1;
alpha = 1 / (n_years * 365.25 * 24);
[fm, x1_1Yr, x2_1Yr] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);

figContour = figure();
plot(tz, hs, '.k');
hold on
plot(x2_1Yr{1}, x1_1Yr{1}, '--r');


n_years = 50;
alpha = 1 / (n_years * 365.25 * 24);
[fm, x1_50Yr, x2_50Yr] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);
figure(figContour);
plot(x2_50Yr{1}, x1_50Yr{1}, '-r');

n_years = 500;
alpha = 1 / (n_years * 365.25 * 24);
[fm, x1_500Yr, x2_500Yr] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);
figure(figContour);
plot(x2_500Yr{1}, x1_500Yr{1}, ':b');

xlabel('Zero-up-crossing period (s)');
ylabel('Significant wave height (m)');


legend({'Observations', '1-year contour', '50-year contour', '500-year contour'}, ...
    'location', 'northwest');
legend box off
box off


% PMS.name = 'Hs-Steepness model';
% PMS.modelType = 'CMA';
% PMS.distributions = {'exponentiated-weibull'; 'lognormal'};
% PMS.isConditionals = {[0 0 0 ]; [1 1]};
% PMS.coeffs = {
%     {pdV_ExpWbl.Alpha pdV_ExpWbl.Beta pdV_ExpWbl.Delta}; 
%     { 
%     @(x1)fMuS(x1);
%     @(x1)fSigmaS(x1)}
%     };
% PMS(1).labels = {'Sign. wave height (m)';
%     'Stepness (-)'};
% PMS.gridCenterPoints = {0.05:0.1:50; 0.001:0.001:0.15};
% 
% n_years = 1;
% alpha = 1 / (n_years * 365.25 * 24);
% [fm, x1_50Yr, x2_50Yr] = computeHdc(PMS, alpha, PMS.gridCenterPoints, 0);
% 
% figure();
% plot(steepness, hs, '.k');
% hold on
% plot(x2_50Yr{1}, x1_50Yr{1}, '-r');
% xlabel('Steepness (-)');
% ylabel('Significant wave height (m)');

%figure(figSteepness)
%x = [0:0.1:max(hs)];
%predictedfSteepnessMedian = fMuTz.a * log(x + fMuTz.b);
%plot(x, predictedfSteepnessMedian)


figure
plot(hs, tz, '.k');
hold on
plot(meanOfHs, medianOfTz, 'or');
predictedTzMedian = exp(fMuTz(x));
plot(x, predictedTzMedian);
xlabel('hs (m)');
ylabel('tz (s)');