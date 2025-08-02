function [out,w] = slcgra(input, calib, weights, method)
% call functions
%       [out,w] = slcgra(input,reference,w,'spsg');
%       [out,w] = slcgra(input,[],w);
% input:    input, rawdata, coil * x * y , sms slice
%           calib, reference data, coil * x * y * z, shifted already
%           weights: output of slcgra also, here should have field
%                ker, [Rx Ry] of the k-space, Rx (odd) * (Ry-1) (even),
%                   supposing acceleration in Ry direction
%                (optional) lambda
%           method, 'spsg' or 'sg', default to be 'spsg'
%  second way:
%           weights, output of the first way, has fields:
%               ker, z, weights
%                   ker, odd in both x, y axes
%                   z, slices in z
%                   weights, coil * (coil*Rx*Ry) * z
%               (optional) lambda
%
if nargin<4
    method = 'spsg';
end
if isfield(weights,'weights')
    out = apply_weights(input,weights);
    w = weights;
else
    switch method
        case 'spsg'
            if isfield(weights,'lambda')
                w   =   weights_spsg(calib, weights.ker,weights.lambda);
            else
                w   =   weights_spsg(calib, weights.ker);
            end
        case 'sg'
            if isfield(weights,'lambda')
                w   =   weights_sg(calib, weights.ker);
            else
                w   =   weights_sg(calib, weights.ker);
            end
    end
    out =   apply_weights(input, w);
end
end
function w = weights_sg(calib, kernel, lambda)

%   Helper function to compute slice grappa kernel weights
%   
%   MChiew
%   Nov 2016

%   calib is [c,kx,ky,z]
%   kernel is (kx, ky), where kx and ky should be odd
%   lambda is an optional regularisation parameter relative to norm(src)

if nargin < 3
    lambda = 0;
end

dims    =   size(calib);
w.ker   =   kernel;
w.z     =   dims(4);

%   Get source and target indices
[src_idx, trg_idx]  =   get_indices(dims(2:3), w.ker);

%   Reshape calib to linearize kx, ky
calib   =   reshape(calib, dims(1), [], dims(4));

%   Select source points on summed calibration data
src =   reshape(sum(calib(:, src_idx, :),3), [], length(trg_idx));

for z = 1:dims(4)
    %   Select target points for each slice
    trg =   calib(:, trg_idx, z);
    
    %   Get (regularised) pseudoinverse
    m   =   src'*pinv(src*src' + norm(src)*lambda*eye(size(src,1)));

    %   Fit weights
    w.weights(:, :, z)  =   trg*m;
end
end
function w = weights_spsg(calib, kernel, lambda)

%   Helper function to compute split slice grappa kernel weights
%   
%   MChiew
%   Nov 2016

%   calib is [c,kx,ky,z]
%   kernel is (kx, ky), where kx and ky should be odd
%   lambda is an optional regularisation parameter relative to norm(src)

if nargin < 3
    lambda = 0;
end

dims    =   size(calib);
w.ker   =   kernel;
w.z     =   dims(4);

%   Get source and target indices
[src_idx, trg_idx]  =   get_indices(dims(2:3), w.ker);

%   Reshape calib to linearize kx, ky
calib   =   reshape(calib, dims(1), [], dims(4));

for z = 1:dims(4)
    %   Select target points for each slice
    trg =   calib(:, trg_idx, z);

    %   Select source points for target slice
    src =   reshape(calib(:, src_idx, z), [], length(trg_idx));

    %   Select all non-target slice locations
    zz  =   setdiff(1:dims(4), z);

    %   Generate zeroed target points
    trg2=   repmat(zeros(size(trg)), 1, length(zz));

    %   Select source points for non-target slice
    src2=   reshape(calib(:, src_idx, zz), [], length(trg_idx)*length(zz));

    %   Get (regularised) pseudoinverse
    m   =   [src src2]'*pinv([src src2]*[src src2]' + norm(src)*lambda*eye(size(src,1)));

    %   Fit weights
    w.weights(:, :, z)  =   [trg trg2]*m;
end
end
function [src, trg] = get_indices(dims, kernel)

%   Helper function to compute source, target and kernel indices
%   
%   MChiew
%   Nov 2016

%   dims is (nx, ny)
%   kernel is (kx, ky), where kx and ky should be odd

pad =   ceil((kernel-1)/2);
ks  =   prod(kernel);% so (ks+1)/2 is the middle ksp target

%   Find calib boundary padding
kx  =   1+pad(1):dims(1)-pad(1);
ky  =   1+pad(2):dims(2)-pad(2);

%   Find relative kernel indices
mask    =   false(dims);
mask(1:kernel(1), 1:kernel(2))  =   true;
k_idx   =   find(mask);
k_idx   =   k_idx - k_idx((ks+1)/2);

%   Find target linear indices
mask    =   false(dims);
mask(kx, ky)    =   true;
trg     =   find(mask);

%   Find source linear indices
src =   bsxfun(@plus, k_idx, trg');
end

function out = apply_weights(input, w)

%   Helper function to apply slice unaliasing grappa weights
%
%   MChiew
%   Nov 2016

%   input is [c,kx,ky,1,t]
%   w is a struct with the kernel and weights ([c,kernel,z])
%   kernel is (kx, ky)

dims    =   size(input);
dims(5) =   size(input,5);
out     =   zeros([dims(1) prod(dims(2:3)) size(w,3) dims(5)]);

%   Pad input boundaries
pad     =   ceil((w.ker-1)/2);
input   =   padarray(input, [0, pad, 0, 0]);

%   Get source and target indices
[src_idx, trg_idx]  =   get_indices(dims(2:3)+2*pad, w.ker);

%   Set periodic boundary condition
input(:, 1:pad(1), :, :, :) =   input(:, end-2*pad(1)+1:end-pad(1), :, :, :);
input(:, end-pad(1)+1:end, :, :, :) =   input(:, pad(1)+1:2*pad(1), :, :, :);
input(:, :, 1:pad(2), :, :) =   input(:, :, end-2*pad(2)+1:end-pad(2), :, :);
input(:, :, end-pad(2)+1:end, :, :) =   input(:, :, pad(2)+1:2*pad(2), :, :);

%   Reshape input to linearize kx, ky
input   =   reshape(input, dims(1), [], 1, dims(5));

for t = 1:dims(5)
    %Select source points for each time point
    src =   reshape(input(:, src_idx, 1, t), [], length(trg_idx));

    %   Apply weights
    for z = 1:w.z
        out(:, :, z, t) =   w.weights(:, :, z)*src;
    end
end

%   Reshape output
out =   reshape(out,[dims(1:3), w.z, dims(5)]);
end