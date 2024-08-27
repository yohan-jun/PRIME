function [ output_args ] = genNii( img, voxel_size, filename )
%GENNII Summary of this function goes here
%   Detailed explanation goes here



save_nii(make_nii(single(img), voxel_size, [], 16), filename)


end

