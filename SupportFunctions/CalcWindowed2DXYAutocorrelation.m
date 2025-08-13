function [acArray, array, varargout] = CalcWindowed2DXYAutocorrelation...
    (array, windowSize, dims, aveDims)
  
  % window must be 2D
  
  % Calcuate the 2D autocorrelation for a 2D array
  % array:        Multidimensional array over which we will calculate
  %               the autocorrelation in windows in 2D across 2 given dimensions,
  %               while averaging over the other dimensions. 
  % windowSize:   [windowSizeX, windowSizeY]  
  % dims:         [ACdimension1, ACdimension2] May be x and y dimensions. 
  % aveDims       [averagingDimension] May be z dimension.
  % 
  % 
  % EXAMPLE
  % [acDeltaXY, array, acDeltaXYMeanMag] = CalcWindowed2DXYAutocorrelation(...
  %   thisDisplacementDS, [2,2] , [2,3], [1]);
  %
  %
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1,2) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>

  % Get dimensions of interest as 1st and 2nd index.
  nDims = ndims(array);
  permuteVec1 = 1:nDims;
  permuteVec1(1) = dims(1);
  permuteVec1(dims(1)) = 1;
  array = permute(array, permuteVec1);
  permuteVec2 = 1:nDims;
  permuteVec2(2) = dims(2);
  permuteVec2(dims(2)) = 2;
  array = permute(array, permuteVec2);
  arrayDims = size(array);
  colonOp = repmat({':'}, 1, nDims);
  nWindowsDim1 = floor(arrayDims(1) / windowSize(1));
  nWindowsDim2 = floor(arrayDims(2) / windowSize(2));
  
  % Crop array if a whole number of windows don't fit.
  if nWindowsDim1 * windowSize(1) ~= arrayDims(1) || nWindowsDim2 * windowSize(2) ~= arrayDims(2)
    array = array(1:nWindowsDim1 * windowSize(1), 1:nWindowsDim2 * windowSize(2), colonOp{3:end});
  end
  
  % Now reshape array. nWindows dimension should always come right after
  % the dimension we are reslicing, so we will not get the right shape we
  % want just yet.
  array = reshape(array, windowSize(1), nWindowsDim1, windowSize(2), nWindowsDim2, arrayDims(3:end));
  % Now we permute.
  array = permute(array, [1 3 2 4:nDims + 2]);
  % Now we have windows of the desired size along dims 1 and 2, then
  % windows indexes in dims 3 and 4, then the rest of the original dims after.
  
  % The AC indices are in dims. aveDims are indices over
  % which the AC will be averaged (most likely 1 or z). 
  % Now, because aveDims is defined wrt the original shape, we need to
  % add 2 dimensions to any dimensions > 2 plus 2 for the additional window
  % dimensions.
  aveDims(aveDims <= 2) = aveDims(aveDims <= 2) + 4;
  acDimSizes = size(array, [1 2]);
  padSize = zeros(1, nDims);
  padSize(1) = acDimSizes(1) - 1;
  padSize(2) = acDimSizes(2) - 1;

  % Compute magnitude for masking.
  if nargout > 1
    % We are asking for mean magnitude inside each window.
    acWinMeanMag = mean(abs(array), [1 2 aveDims]);
    varargout{1} = acWinMeanMag;
  end
  
  % Now compute 2D autocorrelation for each window in the Fourier domain. 
  % Zero-pad the signal for non-circular autocorrelation.
  acArray = abs(...
    fft(fft(...
    padarray(array, padSize, 'post')...
    , [], 1), [], 2)...
    ) .^ 2;
  if ~isempty(aveDims)
    acArray = mean(acArray, aveDims);
  end
  acArray = fftshift(fftshift(ifft(ifft(acArray, [], 1), [], 2), 1), 2);
  
  % Now remove bias by dividing by number of actual elements for each
  % displacement.
  biasDim1 = shiftdim(triang(2 * acDimSizes(1) - 1), 1 - 1) * windowSize(1);
  biasDim2 = shiftdim(triang(2 * acDimSizes(2) - 1), 1 - 2) * windowSize(2);
  biasMat = biasDim1 .* biasDim2;
  acArray = acArray  ./ biasMat;

  % Permute back dimensions.
  if all(dims == [2 3])
    acArray = permute(acArray, [5:nDims + 2, 1:4]);
  end
  
end