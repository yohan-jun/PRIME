function [outputArg1,outputArg2] = write_acq_params(file_name, acq_par)
%WRITE_ACQ_PARAMS Summary of this function goes here
%   Detailed explanation goes here



fid = fopen(file_name, 'w');

fprintf(fid, '%d %d %d %6.5f', acq_par(1,:));
fprintf(fid, '\n%d %d %d %6.5f', acq_par(2,:));

fclose(fid);

end

