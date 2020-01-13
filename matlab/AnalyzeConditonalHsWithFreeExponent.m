DATASET_CHAR = 'D';
BIN_WIDTH = 2;
MIN_DATA_POINTS_IN_BIN = 200;

load datasets-provided-ABCDEF.mat

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

nOfBins = ceil((max(v) / BIN_WIDTH));
pdHs_Exp = ExponentiatedWeibull.empty(nOfBins, 0);
xHat_Exp = cell(nOfBins, 1);
xi = cell(nOfBins, 1);
pstar_i = cell(nOfBins, 1);
xHat_2pWbl = cell(nOfBins, 1);
xHat_Exp = cell(nOfBins, 1);
hsInBin = cell(nOfBins,1);
lowerLimit = nan(nOfBins, 1);
upperLimit = nan(nOfBins, 1);
meanOfHs = nan(nOfBins, 1);
stdOfHs = nan(nOfBins, 1);


fig = figure('position', [100, 100, 1400, 350], 'renderer', 'Painters');
for i = 1:nOfBins
    lowerLimit(i) = (i - 1) * BIN_WIDTH;
    upperLimit(i) = i * BIN_WIDTH;
    vIsInBin = (v > lowerLimit(i)) .* (v < upperLimit(i));
    hsInBin{i} = hs(logical(vIsInBin));
    

    
    if length(hsInBin{i}) >= MIN_DATA_POINTS_IN_BIN
        xi{i} = sort(hsInBin{i});
        n = length(xi{i});
        j = [1:1:n]';
        pi = (j - 0.5) / n;
        pstar_i{i} = computePStar(pi, 1);
        
        pdHs_Exp(i) = ExponentiatedWeibull();
        pdHs_Exp(i).fitDist(hsInBin{i}, 'WLS');
        parmHat = wblfit(hsInBin{i});
        pdHs_2pWbl(i) = makedist('Weibull', parmHat(1), parmHat(2));

        xHat_2pWbl{i} = pdHs_2pWbl(i).icdf(pi);
        xHat_Exp{i} = pdHs_Exp(i).icdf(pi);
        
        % Compute some metrics for analysis
        meanOfHs(i) = mean(hsInBin{i});
        stdOfHs(i) = std(hsInBin{i});
    end
end

meanOfHs = meanOfHs(~isnan(meanOfHs));
stdOfHs = stdOfHs(~isnan(stdOfHs));

for i = 1:length(pdHs_Exp)
        subplot(2, ceil(length(pdHs_Exp) / 2), i);
        plot(log10(xi{i}), pstar_i{i}, 'xk', 'markersize', 4);
        hold on
        plot(log10(xHat_2pWbl{i}), pstar_i{i}, '-k');
        plot(log10(xHat_Exp{i}), pstar_i{i}, '-r');
        xticks = [0.2 0.5 1 2 4 6 8 10 12];
        set(gca, 'xtick', log10(xticks));
        set(gca, 'xticklabel', xticks);
        yticks = [0.05 0.1 0.2 0.5 0.9 0.99 0.9999];
        set(gca, 'ytick', computePStar(yticks, 1));
        set(gca, 'yticklabel', yticks);
        if i > ceil(length(pdHs_Exp) / 2)
            xlabel('Significant wave height (m)');
        end
        if mod(i, ceil(length(pdHs_Exp) / 2)) == 1 
            ylabel('Probability (-)');
        end
        ylim(computePStar([0.05 0.999999], 1))
        grid on
        title([num2str(lowerLimit(i)) ' < v < ' num2str(upperLimit(i)) ' m/s']);
end

alphas = nan(length(pdHs_Exp), 1);
betas = nan(length(pdHs_Exp), 1);
deltas = nan(length(pdHs_Exp), 1);
binCenter = nan(length(pdHs_Exp), 1);
for i = 1:length(pdHs_Exp)
    alphas(i) = pdHs_Exp(i).Alpha;
    betas(i) = pdHs_Exp(i).Beta;
    deltas(i) = pdHs_Exp(i).Delta;
    binCenter(i) = (i - 0.5) * BIN_WIDTH;
end


figure
subplot(3,1,1);
plot(binCenter, alphas, 'ok');
ylabel('\alpha_{hs}');

subplot(3,1,2);
plot(binCenter, betas, 'ok');
ylabel('\beta_{hs}');

subplot(3,1,3);
plot(binCenter, deltas, 'ok');
ylabel('\delta_{hs}');
xlabel('wind speed (m/s)');

figForPaper = figure('position', [100 100 200 180]);
marker_face_color = [0.8 0.8 0.8];
plot(binCenter, deltas, 'ok', 'markerfacecolor', marker_face_color);
hold on
plot([0 30], [5 5], '--k')
ylim([0 25])
xlim([0 30])
box off
xlabel('wind speed (m/s)');
ylabel('\delta_{hs}');
title(['Dataset ' DATASET_CHAR]);

function pstar_i = computePStar(pi, delta)
    pstar_i = log10(-log(1 - pi.^(1 ./ delta)));
end