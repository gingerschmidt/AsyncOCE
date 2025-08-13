function [displacement, displacementPhasor] = CalcDisplacementFromTom(tom, ...
    dopplerFactor, varargin)

  % CALCULATE DISPLACEMENT FROM COMPLEX TOMOGRAM 
  % dopplerFactor:      In µm (wavelengthOCT / (4 * pi * refractiveIdx))
  % bscanGroupSize:     0 = use initial Bscan as reference for everything
  %                     n = reset reference every n frames
  % doSpatialFiltering: true or false, anisotropic complex-averaging to remove noise
  % filterSize:         [windowZ, windowX]
  % 
  % EXAMPLE
  % displacement = displacementFromTom(tom, wavelengthOCT, refractiveIdx,...
  %    'bscanGroupSize', 0, 'doSpatialFiltering', true, 'filterSize', [10,4]);
  %
  %
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1,2) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  % Set default parameters
  bscanGroupSize = 0;
  doSpatialFiltering = 0;
  windowX = 2;
  windowZ = 4;
  % Set optional parameters from varargin
  if(~isempty(varargin))
    for c=1:2:length(varargin)
      switch varargin{c}
        case {'bscanGroupSize'}
          bscanGroupSize = varargin{c+1};
        case {'doSpatialFiltering'}
          doSpatialFiltering = varargin{c+1};
        case {'filterRadii'}
          filterSize = varargin{c+1};
          windowZ = filterSize(1);
          windowX = filterSize(2);
        otherwise
          error(['Invalid optional argument, ', varargin{c}]);
      end
    end
  end
  
  % Calculate displacement relative to first Bscan
  if bscanGroupSize == 0
    displacementPhasor = sum(tom.* conj(tom(:,:,1,:)), 4);
    if doSpatialFiltering
      displacementPhasor = sum(tom.* conj(tom(:,:,1,:)), 4);
      aveKernel = AnisotropicGaussianExp2Diameter([windowX*2+1, windowZ*2+1],windowX*2, windowZ*2);
      displacementPhasor = imfilter(displacementPhasor, aveKernel, 'replicate');
    end
    
  % Calculate displacement from batches of Bscans
  else
    tomSize = size(tom);
    % Set safe group size
    bscanGroupSize = min(bscanGroupSize, tomSize(3));
    tomGroups = reshape(tom(:, :, 1:floor(tomSize(3)/bscanGroupSize)*bscanGroupSize, :),...
      tomSize(1), tomSize(2), bscanGroupSize, [], tomSize(4));
    displacementPhasor = sum(tomGroups(:, :, 2:bscanGroupSize, :, :) .*...
      conj(tomGroups(:, :, 1, :, :)), 5);
    displacementPhasor = reshape(displacementPhasor, tomSize(1), tomSize(2), []);
    if doSpatialFiltering
      aveKernel = AnisotropicGaussianExp2Diameter([windowX*2+1, windowZ*2+1], windowX*2, windowZ*2);
      displacementPhasor = imfilter(displacementPhasor, aveKernel, 'replicate');
    end
  end
  
  % Convert to real units
  displacement = angle(displacementPhasor);
  displacement = displacement * dopplerFactor; % Doppler OCT equation but without dT, in µm
end