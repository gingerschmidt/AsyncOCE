%% ----------------------------------------------------------------------
% To process your own data, you should provide a variable tom (supplied
% here by tom.mat) that contains a complex-valued tomogram with dimensions 
% [Z, X, Bscans, Polarization Channels] AND update other parameters.

load tom.mat
addpath(genpath('SupportFunctions'));

%% ---------------------------------------------------------------------
% Define data parameters
% MODIFY THIS SECTION FOR YOUR OWN DATA
% ----------------------------------------------------------------------

% Define generally constant OCT system parameters
alineRate = 100082;       % In Hz; May need to be calibrated for individual systems.
                          % Precision is critical for out-of-plane correction!
wavelengthOCT = 1.3;      % In µm; wavelength of the OCT source
refractiveIdx = 1.4;      % Refractive index
logLim = [70, 120];       % In dB; limits for tomogram intensity plotting, adjust for your tomogram dynamic range

% Define scanning parameters
nAlinesPerBscan = 1536;   % *BEFORE* any cropping, from original data (important for demodulation)
nBscanStepSize = 2;       % Number of B-scans per y-location
scanWidthX = 12;          % In mm; fast-axis galvo scan range
scanWidthY = 12;          % In mm; slow-axis galvo scan range
noiseFloorROI =...
  {24:32, 500:700};       % Define a region {z0:zEnd, x0:xEnd} in air to calculate the noise floor for intensity masking

% Define elastography parameters
excitationFreq = 1000;    % In Hz; chosen to maximize displacement between B-scans
demodFiltHalfWidth = 12;  % In px; 
rho = 1000;               % In kg/m^3; 

% Surface detection, flattening, and tissue mask parameters
surfaceInterfaceDb = 88;  % In dB; tune this value if surface detection is failing
bigRegionFactor = 11;     % >= 3 for uniformly strong surfaces, lower for surfaces with breaks
maxSurfaceGradient = 30;  % In px/px; used to detect artifactual abrupt changes in tissue height
dcMask = [];              % In your tomogram has a strong DC artifact, provide range so it is removed during surface detection
filterKernel = ...
  ones([31, 15, 1]);      % In px; kernel used for filtering the detected surface. Tune this value if surface detection is failing
filterSizeSurface = ...
  [9 3];                  % In px; kernel used for filtering the detected surface. Tune this value if surface detection is failing
flatSufaceDepth = 100;    % In px; the depth of the flattened tissue surface

% Displacement calculation and visualization parameters
displacementLim = ...
  0.2 * [-1, 1];          % In µm; for displacement visualizations
dispFilterRadii = [4, 0]; % In px; size of Gaussian kernel radii for spatially filtering displacement B-scans
dispMaskFiltSize = ...
  [21, 21, 3];            % In px; size of filter used for masking the displacement

% Demodulation visualization parameters
showBscanDisplacementVideo = true;
showEnFaceDisplacementVideo = true;
nPeriods = 1;             % Number of periods to generate
nTimeStepsPerPeriod = 20; % Number of time steps to generate, per period
generateThisY = 50;       % Generate displacements for this y location (< nY)
  
% Directional Phase Gradient Analysis (DPGA) parameters
gradientDeltaXYDeltaPixel...
  = 1;                    % In px; number of pixels from center over which we compute the phase gradient
zAvgROI = flatSufaceDepth ...
  + (5:30);               % In px; the first few pixels in depth from the surface are limited to Rayleigh surface waves 
airROI =...
  (flatSufaceDepth - 50):...
  (flatSufaceDepth - 10); % In px; air is defined here as the region of the tomogram above the tissue surface after flattening
                          % Used to compute SNR
nBinsForDirFilt = 32;     % Number of directions for directional filtering
corrWindowSizeDeltaXY = ...
  [2, 2];                 % In px; window size for 2D autocorrelation
dims2DAC = [2 3];         % Dimensions over which we will compute the 2D autocorreltion
aveDims2DAC = 1;          % Other dimensions, over which we will be averaging

%% ---------------------------------------------------------------------
% Compute other data parameters
% ----------------------------------------------------------------------

% Tomogram cropping for region of interesting (ROI)
[nZ, nX, nBscans, nChannels] = size(tom);
nY = nBscans / nBscanStepSize;
zROI = 1:nZ;
xROI = 1:nX;
bscanROI = 1:nBscans;

% Crop tomogram
tom = tom(zROI, xROI, bscanROI, :);
% Calculate noise floor
noiseFloorDb = 10 * log10(mean(abs(tom(noiseFloorROI{:}, :, :)) .^ 2, [1 2 3]));

% The demodulation carrier frequency is a function excitation frequency, number
% of samples, and the laser sweep repition rate. See article: 
% Ginger Schmidt, Brett E. Bouma, and Néstor Uribe-Patarroyo, ...
% "Asynchronous, semi-reverberant elastography," Optica 11, 1285-1294 (2024)
demodulationShiftPx = excitationFreq .* nX / alineRate;

% Visualization & graphing settings 
LATEX_DEF = {'Interpreter', 'latex'};
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
diffCMap = colorcet('C3');
cMapDisp = colorcet('D6');
cMapDispMag = viridis;

%% ---------------------------------------------------------------------
% Surface detection & flattening
% ----------------------------------------------------------------------
% Tissue flattening includes surface detection. It also enables us to
% conveniently average across depth to visualize en face surface shear waves. 
tissueFlatteningOptions = struct('surfaceInterfaceDb', surfaceInterfaceDb,...
  'bigRegionFactor', bigRegionFactor, 'dcMask', dcMask,...
  'nBscanStepSize', nBscanStepSize, 'filterKernel', filterKernel,...
  'filterSizeSurface', filterSizeSurface, 'maxSurfaceGradient', maxSurfaceGradient,...
  'desiredTissueSurfaceDepth', flatSufaceDepth, 'figNum', 1);
fprintf('Surface detection...\n');
[tomFlatX, targetSurfaceZPerY, ~, ~] = ...
  DetectTissueSurfaceAndFlatten(tom, tissueFlatteningOptions);
fprintf('done.\n');

%% ---------------------------------------------------------------------
% Compute displacement, unwrap, & surface wave correction
% ----------------------------------------------------------------------
% Compute displacement from tomogram (in um) with Doppler OCT. 
% We assume tissue all has the same refractive index.
dopplerFactor = wavelengthOCT / (4 * pi * refractiveIdx); % To convert from phase to µm
[displacementRaw, displacementPhasorRaw] = ...
  CalcDisplacementFromTom(tomFlatX, dopplerFactor,...
  'bscanGroupSize', nBscanStepSize, 'doSpatialFiltering', true, ...
  'filterRadii', dispFilterRadii);

% Unwrap phase, following Ghiglia DC and Romero LA, Optica Publishing Group (1994)
fprintf('Unwrapping %d Bscans...', size(displacementRaw, 3));
displacementUnwrappedX = (Perform2DDCTPhaseUnwrapping(displacementRaw / dopplerFactor) * dopplerFactor);
fprintf('done.\n');

% Surface wave correction, following Song et. al., J. Biomed. Opt 18, 121505 (2013).
% To avoid propagating displacement to the noise, we only add this corrections
% to regions with sufficient tomogram signal: 6 dB above noise floor and 
% below the tissue surface.
fprintf('Applying surface wave correction...\n');
displacementBscanStepSize = nBscanStepSize - 1;
displacementSurfWaveCorr = displacementUnwrappedX;
displacementMask = 10 * log10(abs(displacementPhasorRaw)) > 10 * log10(mean(10 .^ (noiseFloorDb(:) / 10)));
displacementMask = imclose(displacementMask, strel('cuboid', dispMaskFiltSize));
for thisY = 1:nY
  thisTissueSurface = targetSurfaceZPerY(1, thisY);
  surfaceWaveOffset = (refractiveIdx - 1) *...
    displacementUnwrappedX(thisTissueSurface + 10, :, thisY);
  displacementSurfWaveCorr(thisTissueSurface:end, :, thisY) =...
    displacementUnwrappedX(thisTissueSurface:end, :, thisY) +...
    displacementMask(thisTissueSurface:end, :, thisY) .* surfaceWaveOffset;
end
% Unwrap phase again after surface wave correction
fprintf('Unwrapping %d Bscans after surface wave correction...', size(displacementRaw, 3));
displacementSurfWaveCorr = (Perform2DDCTPhaseUnwrapping(displacementSurfWaveCorr / dopplerFactor) * dopplerFactor);
fprintf('done.\n');

%% ---------------------------------------------------------------------
% In-plane demodulation 
% ----------------------------------------------------------------------
% In the following section, we demodulate the effect of raster-scanning to
% recover the coherent displacement field within B-scans (XZ plane).
% [Refer to sections 2.1.2 and 2.1.4 for further details]. 
fprintf('Demodulating shear wave displacent field... \n');
% Create phase ramp for amplitude demodulation along x, centered on zero. 
xVect = -fix(nX/2):ceil(nX/2) - 1;
% The slope of the phase ramp is f0/fA
ampDemodPhaseSlope = excitationFreq/alineRate;
% Apply frequency domain shift 
displacementDemod = displacementUnwrappedX .* exp(2i * pi * ampDemodPhaseSlope * xVect);
% There are still two copies of the signal (one at 0 and the other at 2*fc)
% Apply a lowpass filter to remove the second copy and any other noise. 
displacementDemodFFT =  fftshift(fft(fftshift(displacementDemod, 2), [], 2), 2);
% Average across y for mean power spectrum after shifting, prior to filtering.
displacementDemodPS_prefilter = mean(abs(displacementDemodFFT) .^ 2, 3);
% Implement filter
demodFilter = blackman(2 * demodFiltHalfWidth + 1);
demodFilter = padarray(demodFilter, nX/2 - demodFiltHalfWidth, 'pre');
demodFilter = padarray(demodFilter, nX/2 - demodFiltHalfWidth - 1, 'post');
demodFilter = demodFilter.';
displacementDemodFFT = displacementDemodFFT .* demodFilter;
displacementDemod = ifftshift(ifft(ifftshift(displacementDemodFFT, 2), [], 2), 2);
% Average across y for mean power spectrum after filtering
displacementDemodPS = mean(abs(displacementDemodFFT) .^ 2, 3);

% ---------------------------------------------------------------------
% Out-of-plane, 3D coherent wave recovery
% ----------------------------------------------------------------------
% The displacement field is now coherent within B-scans. 
% We correct for out-of-plane timing offsets by multiplying with a phase shift
% to account for the time that passes between y-locations. 
% [Refer to sections 2.1.3 and 2.1.4 for further details]. 
% Create complex exponential phase term
yVect = 0:nY-1;
deltaTBetweenYSteps = nBscanStepSize*nAlinesPerBscan/alineRate;
phaseDiffBetweenYSteps = angle(exp(1i*2*pi*(deltaTBetweenYSteps).*yVect*excitationFreq));
% Move the phase differences to the third dimension.
phaseDiffBetweenYSteps = permute(phaseDiffBetweenYSteps, [1 3 2]);
% Apply out-of-plane correction
displacementDemod3D = displacementDemod .* exp(1i .* phaseDiffBetweenYSteps );
% Displacement is phase stable along y now, so we can apply spatial filtering. 
displacementDemod3DFFT = fftshift(fft(fftshift(displacementDemod3D, 3), [], 3), 3);
displacementDemod3DFFT_unfiltered = displacementDemod3DFFT;
displacementYFilter = blackman(nY);
displacementYFilter = permute(displacementYFilter, [3 2 1]);
displacementDemod3DFFT = displacementDemod3DFFT.* displacementYFilter;
displacementDemod3D = ifftshift(ifft(ifftshift(displacementDemod3DFFT, 3), [], 3), 3);
% Done with downshifting and filtering
fprintf('done.\n');

% Visualize demodulation power spectra
% Power spectrum for in plane demodulation, before and after filtering.
figure(2); subplot(2,1,1)
psFreqArray = -round(demodulationShiftPx) * 3:round(demodulationShiftPx) * 3;
imagesc(psFreqArray, 1:nZ, 10 * log10(displacementDemodPS_prefilter(:, nX/2 + 1 + psFreqArray) .^ 2)), colormap(gca, colorcet('L20')), colorbar
hAx = gca;
psLims = hAx.CLim;
title({'In-plane, B-scan power spectrum taken along x and averaged across y', ...
  'After shifting and before filtering'}), ylabel('Depth');
subplot(2,1,2)
imagesc(psFreqArray, 1:nZ, 10 * log10(displacementDemodPS(:, nX/2 + 1 + psFreqArray) .^ 2), psLims), colormap(gca, colorcet('L20')), colorbar
title('After filtering'), xlabel('log10 Power spectrum zoomed-in'), ylabel('Depth');
% Power spectrum for before and after out-of-plane filtering along y.
figure(3); subplot(2,1,1)
psFreqArray = -round(demodulationShiftPx) * 5:round(demodulationShiftPx) * 5;
imagesc(10 * log10(abs(squeeze(mean(displacementDemod3DFFT_unfiltered,1))).^ 2)), colormap(gca, colorcet('L20')), colorbar
hAx = gca;
psLims = hAx.CLim;
title({'En face power spectrum taken along y:', 'Before filtering'}), ylabel('Depth');
subplot(2,1,2)
imagesc(10 * log10(abs(squeeze(mean(displacementDemod3DFFT,1))).^ 2), psLims), colormap(gca, colorcet('L20')), colorbar
title('After filtering'), xlabel('log10 Power spectrum zoomed-in'), ylabel('Depth');

%% ---------------------------------------------------------------------
% Generate displacement at *any* time step from just 2 Bscans at each y
% ----------------------------------------------------------------------
% Compute phase array
thisPhaseVec = linspace(0, 2 * pi * nPeriods, nPeriods * nTimeStepsPerPeriod);
% Plot B-scan (XZ plane) displacement video
if showBscanDisplacementVideo
  for thisPhase = thisPhaseVec
    figure(4);
    imagesc(squeeze(real(displacementDemod(:, :, generateThisY)...
      .* exp(-1i .* thisPhase))), displacementLim), colorbar, colorbarlbl('[\mum]'), colormap(cMapDisp)
    title('Demodulated B-scan (XZ plane) displacement'), xlabel('x'), ylabel('z')
    drawnow
  end
end
% Plot en face (XY plane) displacement video
if showEnFaceDisplacementVideo
  for thisPhase = thisPhaseVec
    figure(5);
    imagescnan(squeeze(real(squeeze(mean(displacementDemod3D(zAvgROI,:,:),1)).' .*...
      exp(-1i .* thisPhase))), displacementLim), colorbar, colorbarlbl('[\mum]'), colormap(cMapDisp)
    title('Demodulated en face displacement'), xlabel('x'), ylabel('y')
    drawnow
  end
end

%% -------------------------------------------------------------
% Directional filtering
% -------------------------------------------------------------
% First, downscale in X to match Y.
nYDS = nY;
nXDS = nYDS;
nZDS = numel(zAvgROI);

% Downsample
displacementDS = ipermute(imresize(permute(displacementDemod3D, ...
  [2 3 1]), [nXDS nYDS]), [2 3 1]);

% Calculate NoiseFloor
displacementNoiseFloor = mean(abs(displacementDS(airROI, :, :)), [1 2 3]);

% Compute 2D FFT across x and y dimensions, which should also be the same length now. 
displacementDSFTXY = fftshift(fftshift(fft(fft(ifftshift(ifftshift(displacementDS(zAvgROI, :, :), 2), 3), [], 2), [], 3), 2), 3);

% Create directional filters
[dirFiltBank] = CreateDirectionalFilters(nBinsForDirFilt, nXDS, nYDS, 6);

% Apply directional filters to entire displacement field
displacementDSFTDirFilt = displacementDSFTXY .* dirFiltBank;
displacementDSDirFilt = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(displacementDSFTDirFilt, 2), 3), [], 2), [], 3), 2), 3);

% Visualize directional filtering results
figure(7); clf;
tiledlayout(3,2,'Padding', 'none', 'TileSpacing', 'compact');
for thisDir = 5:5:nBinsForDirFilt
  figure(7), nexttile, imagescnan(squeeze(real(mean(displacementDSDirFilt(:, :, :, thisDir), 1) .*...
      exp(-1i .* thisPhase))).', displacementLim * 0.2), colorbar, colormap(cMapDisp)
  title(['Displacement ' num2str(thisDir)]); xlabel('x [px]'); ylabel('y [px]');
  drawnow
end
%% -------------------------------------------------------------
% AC phase gradient
% ------------------------------------------------------------
xAxisSpacing = scanWidthX / nXDS; % in mm/px
yAxisSpacing = scanWidthY / nYDS; % in mm/px
nWindows2DX = floor(nXDS / corrWindowSizeDeltaXY(1));
nWindows2DY = floor(nYDS / corrWindowSizeDeltaXY(2));
% Calculate autocorrelation
acDeltaXYPhaseChangeArray = zeros([nWindows2DX, nWindows2DY, nBinsForDirFilt], 'single');
acDeltaXYMeanMagArray = zeros([nWindows2DX, nWindows2DY, nBinsForDirFilt], 'single');
% Loop through each directionally filtered displacement field
for thisDir = 1:nBinsForDirFilt
  thisDisplacementDS = displacementDSDirFilt(:, :, :, thisDir);
  acDeltaXYPhaseChange = zeros([nWindows2DX, nWindows2DY], 'single');
  [acDeltaXY, array, acDeltaXYMeanMag] = CalcWindowed2DXYAutocorrelation(...
    thisDisplacementDS, corrWindowSizeDeltaXY, dims2DAC, aveDims2DAC);
  % Normalize ACs
  acDeltaXY = acDeltaXY ./ acDeltaXY(:, corrWindowSizeDeltaXY(1), corrWindowSizeDeltaXY(2), :, :);
  % Compute phase gradient within each window
  for thisXWin = 1:nWindows2DX
    for thisYWin = 1:nWindows2DY
      thisAC = shiftdim(acDeltaXY(1, :, :, thisXWin, thisYWin), 1);
      thisACPhase = angle(thisAC);
      thisACPhase = Perform2DDCTPhaseUnwrapping(thisACPhase);
      % All Xs and Ys here will be flipped from the expected order because the
      % 1st index is X and second index is Y which is the opposite of how gradient works.
      [acPhaseGradY, acPhaseGradX] = gradient(thisACPhase, yAxisSpacing, xAxisSpacing);
      deltaXAxis = -corrWindowSizeDeltaXY(1) + 1:corrWindowSizeDeltaXY(1) - 1;
      deltaYAxis = -corrWindowSizeDeltaXY(2) + 1:corrWindowSizeDeltaXY(2) - 1;
      acPhaseGradMag = sqrt(acPhaseGradX(corrWindowSizeDeltaXY(1), corrWindowSizeDeltaXY(2)) .^ 2 +...
        acPhaseGradY(corrWindowSizeDeltaXY(1), corrWindowSizeDeltaXY(2)) .^ 2);
      acDeltaXYPhaseChange(thisXWin, thisYWin) = acPhaseGradMag;
      acDeltaXYPhaseChangeArray(:, :, thisDir) = acDeltaXYPhaseChange;
      acDeltaXYMeanMagArray(:, :, thisDir) = acDeltaXYMeanMag;
    end
  end
end

%% ---------------------------------------------------------------------
% AsyncOCE and DPGA RESULTS
% ----------------------------------------------------------------------
% Mask by magnitude relative to displacement noise floor SNR.
acDeltaXYPhaseChangeArrayMasked = acDeltaXYPhaseChangeArray;
acDeltaXYPhaseChangeArrayMasked(acDeltaXYMeanMagArray < displacementNoiseFloor) = nan;
% Average across all directions, assuming a laterally isotropic sample
wavespeedAvgEnFace = 2 * pi ./ (mean(acDeltaXYPhaseChangeArrayMasked, 3, 'omitnan').') / 1000 * excitationFreq;
shearModAvgEnFace = wavespeedAvgEnFace .* wavespeedAvgEnFace * rho / 10 ^ 3;
% Compute average shear wave number, wave speed, and shear modulus 
% assuming a homogenous sample
kMean = mean(acDeltaXYPhaseChangeArrayMasked, [1,2,3], 'omitnan');
kStd = std(acDeltaXYPhaseChangeArrayMasked, 0, [1,2,3], 'omitnan');
cMean = (2 * pi / kMean) / 1000 * excitationFreq;

%% Plot results
% Cropping out 1 px along each edge due to filtering artifacts
figure(8); imagescnan(mean(acDeltaXYPhaseChangeArrayMasked(2:end - 1,2:end - 1, :),3,'omitnan').', [0 5]), ...
    colorbar, colorbarlbl('[1/mm]'), colormap(viridis), title('wavenumber k');
       xlabel('x'), ylabel('y');
figure(9); imagescnan(wavespeedAvgEnFace(2:end - 1,2:end - 1), [0 10]), ...
    colorbar,colorbarlbl('[m/s]'), colormap(colorcet('L20')), title('wavespeed c');
       xlabel('x'), ylabel('y')
figure(10); imagescnan(shearModAvgEnFace(2:end - 1,2:end - 1), [0 30]), ...
    colorbar, colorbarlbl('[kPa]'), colormap(colorcet('L6')), title('shear modulus');
       xlabel('x'), ylabel('y')
