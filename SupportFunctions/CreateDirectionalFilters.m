function [dirFiltBank] = CreateDirectionalFilters(nBins, nX, nY, figNum)

  % Make pie slice filters based on arctan angle. 
  dirFiltAngles = linspace(0, 2 * pi - 2 * pi / nBins, nBins) - pi;
  [dirFiltX, dirFiltY] = meshgrid(-fix(nX / 2):ceil(nX / 2) - 1, -fix(nY / 2):ceil(nY / 2) - 1);
  [dirFiltTheta, ~] = cart2pol(dirFiltX, dirFiltY);
  dirFiltTheta = angle(exp(1i * (dirFiltTheta + pi)));
  dirFiltBank = zeros(nX, nY, nBins, 'single');
  
  % Visualize filters and resulting displacement fields.
  figure(figNum); clf;
  for thisBin = 1:nBins
    dirFiltAngleHalfWidth = 2 * pi / (nBins);
    if thisBin < floor(nBins / 2)
      dirFiltBank(:, :, thisBin) = exp(-(angle(exp(1i * (dirFiltTheta + pi))) - (dirFiltAngles(thisBin) + pi)) .^ 2 / dirFiltAngleHalfWidth .^ 2);
    else
      dirFiltBank(:, :, thisBin) = exp(-(dirFiltTheta - dirFiltAngles(thisBin)) .^ 2 / dirFiltAngleHalfWidth .^ 2);
    end
    if mod(thisBin, 5) == 0
      figure(figNum), nexttile, imagesc(dirFiltBank(:,:,thisBin)), colormap(gca, colorcet('L20'))
      title(['Directional Filter ' num2str(thisBin)]); xlabel('kx'); ylabel('ky');
    end
  end
  
  % Delete DC component
  dirFiltBank(fix(nY / 2) + 1, fix(nX / 2) + 1, :) = 0;
  
  % Move directions to 4th index
  dirFiltBank = permute(dirFiltBank, [4 2 1 3]);

end

