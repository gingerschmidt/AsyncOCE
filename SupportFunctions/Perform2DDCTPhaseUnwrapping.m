function [phiUnwrapped, phi] = Perform2DDCTPhaseUnwrapping(psi, weights)
  %Perform2DDCTPhaseUnwrapping 2D phase unwrapping following:
  %   Ghiglia DC, Romero LA. Robust two-dimensional weighted and unweighted
  %   phase unwrapping that uses fast transforms and iterative methods. J Opt
  %   Soc Am A, JOSAA. Optica Publishing Group; 1994 Jan 1;11(1):107–117.
  
  %   Uses the discrete cosine transform to find the LSE solution in 2D
  
  if (nargin < 2) && isreal(psi) % No weights
    % Compute the wrapped wrapped-phase differences along 1st index
    deltaZ = wrapToPi(diff(psi, 1, 1));
    % Set boundary conditions
    deltaZ = padarray(deltaZ, [1 0], 'both');
    % Compute the wrapped wrapped-phase differences along 2nd index
    deltaX = wrapToPi(diff(psi, 1, 2));
    % Set boundary conditions
    deltaX = padarray(deltaX, [0 1], 'both');
    
    % Calculate rho, essentially the known side of the Poisson equation (the
    % second derivatives of psi)
    rho = diff(deltaZ, 1, 1) + diff(deltaX, 1, 2);
    
    % Get the unwrapped phase by solving the Poisson equation
    phiUnwrapped = SolvePoissonWithDCT(rho);
  else
    % Perform preconditioned conjugate gradient iterated method for weighted LSE
    % phase unwrapping
    if ~isreal(psi)
      weights = abs(psi);
      psi = angle(psi);
    else
      if (~all(size(weights, [1 2]) == size(psi, [1 2])))
        error('psi and weight need to have the same 1st and 2nd dimensions');
      end
    end
    % First we need to reduce weight to (nZ-1)*(nX-1), according to the paper,
    % let's take the minimum between adjacent values
    weightsDZ = min(weights(1:end - 1, :, :), weights(2:end, :, :));
    weightsDX = min(weights(:, 1:end - 1, :), weights(:, 2:end, :));
    % Compute the Weighted wrapped wrapped-phase differences along 1st index
    deltaZ = wrapToPi(diff(psi, 1, 1)) .* weightsDZ .^ 2;
    % Compute the Weighted wrapped wrapped-phase differences along 2nd index
    deltaX = wrapToPi(diff(psi, 1, 2)) .* weightsDX .^ 2;
    % Set boundary conditions
    deltaZ = padarray(deltaZ, [1 0], 'both');
    % Set boundary conditions
    deltaX = padarray(deltaX, [0 1], 'both');
    % Calculate rho0, essentially the known side of the Poisson equation (the
    % second derivatives of psi), first iteration
    rhoK = diff(deltaZ, 1, 1) + diff(deltaX, 1, 2);
    % And its norm in each plane
    rho0Norm = sqrt(sum(rhoK .^ 2, [1 2]));    
    
    % Now iterate
    tol = 1e-8;
    k = 0;
    doIter = true;
    phiUnwrapped = zeros(size(psi), 'like', psi);
    while doIter
      zK = SolvePoissonWithDCT(rhoK);
      k = k + 1;
      if k == 1
        % First iteration, direct solution
        pK = zK;
        fprintf('Iteration (current tolerance, %.2g target): ', tol)
      else
        % Update vectors
        betaK = sum(rhoK .* zK, [1 2]) ./ sum(rhoKPrev .* zKPrev, [1 2]);
        pK = zK + betaK .* pK;
      end
      % Save current solution for next iteration
      rhoKPrev = rhoK;
      zKPrev = zK;
      
      % Perform step 6, applying Q operator to rho
      deltaZ = (diff(pK, 1, 1)) .* weightsDZ .^ 2;
      % Compute the Weighted wrapped wrapped-phase differences along 2nd index
      deltaX = (diff(pK, 1, 2)) .* weightsDX .^ 2;
      % Set boundary conditions
      deltaZ = padarray(deltaZ, [1 0], 'both');
      % Set boundary conditions
      deltaX = padarray(deltaX, [0 1], 'both');
      % Calculate rho0, essentially the known side of the Poisson equation (the
      % second derivatives of psi), first iteration
      pKQ = diff(deltaZ, 1, 1) + diff(deltaX, 1, 2);
      % Update alphaK, which is a scalar per 2D plane
      alphaK = sum(rhoK .* zK, [1 2]) ./ sum(pK .* pKQ, [1 2]);
      phiUnwrapped = phiUnwrapped + alphaK .* pK;
      rhoK = rhoK - alphaK .* pKQ;  
      
      % Check for termination condition. We use the mean norm instead of all the
      % norms to avoid getting stuck in the iterations by a single problematic
      % frame
      tolCurr = mean(sqrt(sum(rhoK .^ 2, [1 2])) ./ rho0Norm);
      if (k >= prod(size(psi, [1 2]))) || (tolCurr < tol)
        doIter = false;
      else
        fprintf(' %d (%.2g),', k, tolCurr)
      end
    end
    fprintf(' done.\n,')
  end
end

function phi = SolvePoissonWithDCT(rho)
  % Solve the Poisson equation using dct
  % Apply 2D DCT
  rhoDCTZX = dct(dct(rho, [], 1), [], 2);
  % Calculate vectors and factors needed for the solution. We want these to be
  % in double precision
  [nZ, nX] = size(rho, [1 2]);
  [xMat, zMat] = meshgrid(0:nX - 1, 0:nZ - 1);
  % Now the solution is the 2D DCT multiplied by the corresponding factors
  phiDCTZX = rhoDCTZX ./ (2 .* (cos(pi .* xMat ./ nX) + cos(pi .* zMat ./ nZ) - 2));
  % We preserve the bias (mean) of the wrapped phase by making the undefined 0,0
  % term equal to rhoDCTZX(0,0) (instead of 0 bias like the old function).
  phiDCTZX(1, 1, :) = rhoDCTZX(1, 1, :);
  % Now apply inverse DCT to get the result
  phi = idct(idct(phiDCTZX, [], 1), [], 2);
  
end