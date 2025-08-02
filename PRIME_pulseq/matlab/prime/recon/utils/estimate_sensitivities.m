function [sensitivities] = estimate_sensitivities(data, method, b,cnum)
% function [sensitivities] = estimate_sensitivities(data, method, b,cnum)
%
% Adaptive combination weights estimation:
%
% Input:
%   data: [Nx,Ny,Nc]
%   2D complex images [Nx,Ny] from Nc coils
%   assumed to be prewhitened
%
%   method: integer
%      0: proposed method (default)
%      1: Walsh et al., MRM 2000;43(5)
%
%   b: integer
%      box size=2*b+1
%      default = 3 i.e. box size = 7
%
%   cnum: integer
%      Reference coil number for Walsh
%      default = 1
%
% Output:
%   sensitivities: [Nx,Ny,Nc]
%   Estimated relative coil sensitivities
%
[Nx,Ny,Nc] = size(data);

if (nargin<2)
    method = 0;
    b = 3;
    cnum = 1;
elseif (nargin < 3)
    b = 3;
    cnum = 1;
elseif nargin<4
    cnum = 1;
end

% Initialize the result
sensitivities = zeros(Nx,Ny,Nc);

% Loop over the pixels in the image
% We ignore the edges here for simplicity
% A full implementaion would handle the edges properly.
for x = 1:Nx
  for y = 1:Ny

    % Get the data in a box around the current point
    % Use a 7x7 box.  Size is arbitrary, i.e. coil smoothness assumption
    xmin = max(1,x-b); xmax = min(Nx,x+b);
    ymin = max(1,y-b); ymax = min(Ny,y+b);
    bx = xmax-xmin+1; by = ymax-ymin+1;
    D = reshape(data(xmin:xmax,ymin:ymax,:),[bx*by,Nc]);

    % Compute the sensitivity
    if (method == 0)
      % Proposed method

      % The SVD way
      % [U,S,V]   = svd(D);
      % sorted from biggest to smallest
      % u1 = U(:,1); % first left singular vector (Bx*By,1)
      % v1 = V(:,1); % first right singular vector (Nc,1)
      % s1 = S(1,1); % first singular value
      % the above is slow, and since we only care about the first
      % singular vector, we can use the power method to solve
      
      % Power method
      % initialize to the mean of the data
      v1 = transpose(mean(D,1)); 
      v1 = v1/norm(v1);
      % 3 iterations
      for iter = 1:3
        v1 = D'*D*v1; 
        v1 = v1/norm(v1);
      end
      % compute u1 and s1
      u1 = D*v1;
      s1 = norm(u1);
      u1 = u1/s1;
    
      % Compute the phase of the average of u1
      theta = angle(mean(u1));

      % Normalized coil sensitivies (combination weights)
      sensitivities(x,y,:)  = exp(1i*theta)*v1';
      
    else
      % Original Walsh et al. method

      % Compute the covariance over the box (an Nc by Nc matrix)
      R = D'*D;
      % The eigenvalue way
      % [V,D]   = eig(R);
      % sorted from smallest to biggest
      % v1 = V(:,end); % biggest eigenvector (Nc,1)
      % the above is slow, and since we only care about the first
      % eigenvector, we can use the power method to solve
      
      % Power method
      % initialize to the mean of the covariance
      v1 = transpose(mean(R,1)); 
      v1 = v1/norm(v1);
      % 3 iterations
      for iter = 1:3
        v1 = R*v1; 
        v1 = v1/norm(v1);
      end
      
      % Compute the phase of the reference coil at this location
      reference_coil = cnum;%1;  % arbitrary
      theta = angle(v1(reference_coil));
      
      % Normalized relative coil sensitivies
      sensitivities(x,y,:)  = exp(1i*theta)*v1';
      
    end

  end % Loop over y
end % Loop over x

end
