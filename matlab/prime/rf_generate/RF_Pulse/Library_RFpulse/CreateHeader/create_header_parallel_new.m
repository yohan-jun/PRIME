%WARNING: need to make sure that gradient and rf go to zero at begin and
%end of the pulse, e.g. for spiral need to have rewinder

function create_header_parallel_new(g,b_est,oversample)

DimRF = 2;
DimGradient = 3;


G_length = length(g);
B_length = G_length*oversample;
n_coils = length(b_est)/B_length;

for count = 1:n_coils
    rf(:,count) = b_est(1+(count-1)*B_length:count*B_length);
end

max_rf = max(abs(rf),[],1);

angle_rfs = angle(rf);
index = angle_rfs<0;
angle_rfs(index) = angle_rfs(index)+2*pi;
 
abs_rfs = abs(rf);

%fid = fopen('pTx.float', 'wb');
fid = fopen('pTx.float', 'w');

fprintf(fid, '#%%Header.ArbitraryRFpulse \n');
fprintf(fid,'[pTXPulse]\n');
fprintf(fid,'NUsedChannels = %g\n', n_coils);
fprintf(fid,'DimRF = %g\n', DimRF);
fprintf(fid,'DimGradient = %g\n', DimGradient);
fprintf(fid,'\n\n');

for ChnCount = 1:n_coils
    fprintf(fid,'# ----------\n');
    fprintf(fid,'[pTXPulse_ch%g]\n', ChnCount-1);
    fprintf(fid,'NRFPoints     =  %g\n', B_length);
    fprintf(fid,'NGradientPoints     =  %g\n\n', G_length);
    fprintf(fid,'MaxRF     =  %g\n\n', max_rf(ChnCount));
    for RFcount = 1:B_length
        fprintf(fid,'RF[%g] = %g     %g\n', RFcount, abs_rfs( RFcount, ChnCount), angle_rfs( RFcount, ChnCount));
    end
    fprintf(fid,'\n');
    for Gcount = 1:G_length
        fprintf(fid,'G[%g] = %g     %g     %g\n', Gcount, g(Gcount,1), g(Gcount,2), g(Gcount,3));
    end
    
    fprintf(fid,'\n\n');
end

fclose(fid);