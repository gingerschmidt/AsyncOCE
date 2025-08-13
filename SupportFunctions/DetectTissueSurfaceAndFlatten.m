function [tomFlatX, targetSurfaceZPerY, meanSurfaceZPerY, tissueMask, surfaceZ] =...
    DetectTissueSurfaceAndFlatten(tom, options)
  
  % Detects tissue surface and flattens to average detected surface depth.
  % tom:  dimensions are [Z, X, Bscans, Polarization Channels]
  % See default params below for additional arguments to tune if surface
  % detection is fail. 
  %
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1)
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  
  if nargin == 0
    error('At the least one input is required');
  elseif nargin > 1
    % If options are provided, unpack them
    StructToVars(options);
  end
  
  %% Default parms
  if ~exist('figNum', 'var') %#ok<*NODEF>
    figNum = [];  % If empty, will not specify figure number
  end
  if ~exist('filterKernel', 'var') || isempty(filterKernel)
    filterKernel = ones([20, 30, 1]);  % Tune this value if surface detection is failing
  end
  if ~exist('filterSizeSurface', 'var') || isempty(filterSizeSurface)
    filterSizeSurface = [3 3];  % Tune this value if surface detection is failing
  end
  if ~exist('surfaceInterfaceDb', 'var') || isempty(surfaceInterfaceDb)
    surfaceInterfaceDb = 89;  % Tune this value if surface detection is failing
  end
  if ~exist('noiseFloorDb', 'var') || isempty(noiseFloorDb)
    % This estimation assumes there are some depths domintated by noise. If
    % this estimation is not accurate, provide value for all channels
    noiseFloorDb = min(mean(abs(tom(:, :, 1, :)) .^ 2, 2), [], 1);
  end
  if ~exist('bigRegionFactor', 'var') || isempty(bigRegionFactor)
    bigRegionFactor = 7;  % Tune this value if surface detection is failing
  end
  if ~exist('desiredTissueSurfaceDepth', 'var')
    desiredTissueSurfaceDepth = 200;  % If empty, will be set to average tissue surface depth at each Y (so only X flattening)
  end
  if ~exist('dcMask', 'var')
    dcMask = [];
  end
  if ~exist('maxSurfaceGradient', 'var')
    maxSurfaceGradient = 15;
  end
  if ~exist('nBscanStepSize', 'var') || isempty(nBscanStepSize)
    nBscanStepSize = 2;  % Only used if desiredTissueSurfaceDepth is empty
  end
  
  if ~exist('surfaceZ', 'var')  || isempty(surfaceZ)
    % Find top surface
    tomIntMaskSurface = sum(abs(tom) .^ 2, 4); % Add both polarization channels
    % Remove DC artifact
    tomIntMaskSurface(dcMask, :, :) = sum(10 .^ (noiseFloorDb / 10));
    tissueMask = CreateTissueMask(...
      tomIntMaskSurface, filterKernel, surfaceInterfaceDb, bigRegionFactor);
    [~, topZ] = max(tissueMask, [], 1);
    topZ = squeeze(topZ);

    % Check detected raw surface is meaningful
    topZ(topZ <= 1) = nan;
    
    % Plot surface detection results
    surfaceLim = median(topZ, 'all', 'omitnan') + [-2 2] .* std(topZ, 0, 'all', 'omitnan');
    if isempty(figNum)
      figure();
    else
      figure(figNum);
    end
    subplot(1, 3, 1), imagescnan(topZ.', surfaceLim), colorbar
    title('Detected surface', 'Raw')
    [topZDZ, topZDX] = gradient(topZ);
    topZDZDX = sqrt(topZDZ .^ 2 + topZDX .^ 2);
    surfaceZ = double(topZ);
    surfaceZ(topZDZDX > maxSurfaceGradient) = nan;
    surfaceZ = medfilt2(surfaceZ, filterSizeSurface, 'symmetric');
    subplot(1, 3, 2), imagescnan(surfaceZ.', surfaceLim), colorbar
    title('Detected surface', 'Cleaned up')
    surfaceZ = inpaint_nans(surfaceZ, 2);
    surfaceZ = ceil(surfaceZ);
    subplot(1, 3, 3), imagescnan(surfaceZ.', surfaceLim), colorbar
    title('Detected surface', 'Interpolated')
    drawnow
  else
    fprintf('Using provided surfaceZ for flattening\n')
  end
  % Use surface detection & masking results to flatten tissue
  % Pad array to avoid tissue wrapping around
  padSize = round(size(tom, 1) / 4);
  tomFlatX = padarray(tom, [padSize, 0], 'post');
  % Remove DC artifact
  [nZ, nX, nY, nChs] = size(tom, [1 2 3 4]);
  nY = nY / nBscanStepSize;
  meanSurfaceZPerY = squeeze(round(mean(reshape(surfaceZ, nX, nBscanStepSize, nY), 2)));
  for thisCh = 1:nChs
    tomFlatX(dcMask, :, :, thisCh) = sqrt(10 .^ (noiseFloorDb(:, :, :, thisCh) / 10));
  end
  if isempty(desiredTissueSurfaceDepth)
    targetSurfaceZPerY = round(mean(meanSurfaceZPerY, 1));
  else
    targetSurfaceZPerY = desiredTissueSurfaceDepth * ones(1, nY);
  end
  for thisY = 1:nY
    for thisX = 1:nX
      tomFlatX(:, thisX, nBscanStepSize * (thisY - 1) + (1:nBscanStepSize), :) =...
        circshift(...
        tomFlatX(:, thisX, nBscanStepSize * (thisY - 1) + (1:nBscanStepSize), :),...
        targetSurfaceZPerY(:, thisY) - meanSurfaceZPerY(thisX, thisY), 1);
    end
  end
  % Remove padding
  tomFlatX = tomFlatX(1:nZ, :, :, :);
end