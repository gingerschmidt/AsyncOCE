function [ tissueMask ] = CreateTissueMask(tomInt, filterKernel,...
    surfaceInterfaceDb, bigRegionFactor)

  % CREATE TISSUE MASK BASED ON INTENSITY
  % tomInt:             Intensity tomogram in liner scale (abs(tom).^2)
  % filterKernal:       [Z, X]
  % surfaceInterfaceDb: Tissue intensity at the surface
  % bigRegionFactor:    >= 3 for uniformly strong surfaces, lower for surfaces with breaks
  % 
  % 
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  
  if isempty(bigRegionFactor)
    bigRegionFactor = 5;
  end
  
  tissueMask = tomInt > 10 ^ (surfaceInterfaceDb / 10);
  tissueMask = imfilter(uint8(tissueMask), filterKernel, 'replicate');
  tissueMask = tissueMask > sum(filterKernel(:)) * 0.2;

  for k = 1:size(tomInt, 3)
    thisMask = tissueMask(:, :, k);
    % Fill holes
    thisMask = imfill(thisMask, 'holes');
    % Find connected regions
    groups = bwconncomp(thisMask, 4);
    % Get number of pixels for each region
    numPixels = cellfun(@numel,groups.PixelIdxList);
    % Sort by size
    [numPixelsSorted, sortIdx] = sort(numPixels, 'descend');
    % Identify really big regions, the ones we want
    bigRegions = numPixelsSorted > bigRegionFactor * size(thisMask, 2); 
    
    % Let's analyze the azimuthal extent of each one
    for thisRegionIdx = 1:numel(numPixelsSorted)
      thisRegion = zeros(size(thisMask), 'logical');
      thisRegion(groups.PixelIdxList{sortIdx(thisRegionIdx)}) = 1;
      % Calc mean thickness and width
      meanThickness = mean(sum(thisRegion, 1), 2);
      regionWidth = sum(sum(thisRegion, 1) > 0.8 * meanThickness, 2);
    end
    validRegions = bigRegions;
    nValidRegions = sum(validRegions);
    if nValidRegions == 0
      % Then get widest region
      nValidRegions = 1;
      [~, validRegionsIdx] = max(regionWidth);
    else
      validRegionsIdx = find(validRegions);
    end
    
    % Create mask based on chosen regions
    thisMask(:) = false;
    for thisRegionIdx = 1:nValidRegions
      thisMask(groups.PixelIdxList{sortIdx(validRegionsIdx(thisRegionIdx))}) = 1;
    end
    
    % Now connect them
    thisMask = imfill(thisMask, 'holes');
    
    % Done!
    tissueMask(:, :, k) = thisMask;
  end
end
