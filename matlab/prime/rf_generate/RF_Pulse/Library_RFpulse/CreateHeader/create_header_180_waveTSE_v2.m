function create_header_180_diffusion(rf,g,factor1,factor2,VoltageFactor,Slice_Seperation,Asymetry,RewindGzMoment)
%g is a 3 x column matrix (not mult by 10 yet)
%rf is a complex column vector

max_x = max(g(:,1));
max_y = max(g(:,2));
max_z = max(g(:,3));
max_rf = max(abs(rf(:)));

min_x = min(g(:,1));
min_y = min(g(:,2));
min_z = min(g(:,3));
min_rf = min(rf);

if abs(min_x)> abs(max_x)
    max_x = abs(min_x);
end
if max_x == 0
        max_x = 0.1;
end
if abs(min_y)> abs(max_y)
    max_y = abs(min_y);
end
if max_y == 0
        max_y = 0.1;
end
if abs(min_z)> abs(max_z)
    max_z = abs(min_z);    
end
if max_z == 0
        max_z = 0.1;
end
if abs(min_rf)> abs(max_rf)
    max_rf = abs(min_rf);    
end
if max_rf == 0
    max_rf = 0.1;
end


fidx = fopen('gx1801.float', 'wb');
fidy = fopen('gy1801.float', 'wb');
fidz = fopen('gz1801.float', 'wb');
fidx2 = fopen('gx1802.float', 'wb');
fidy2 = fopen('gy1802.float', 'wb');
fidz2 = fopen('gz1802.float', 'wb');

fwrite(fidx, max_x*10, 'float32');
fwrite(fidx, length(g(:,1)), 'float32');
fwrite(fidx, g(:,1)*10, 'float32');
fwrite(fidx2, max_x*10, 'float32');
fwrite(fidx2, length(g(:,1)), 'float32');
fwrite(fidx2, g(:,1)*10, 'float32');

fwrite(fidy, max_y*10, 'float32');
fwrite(fidy, length(g(:,2)), 'float32');
fwrite(fidy, g(:,2)*10, 'float32');
fwrite(fidy2, max_y*10, 'float32');
fwrite(fidy2, length(g(:,2)), 'float32');
fwrite(fidy2, g(:,2)*10, 'float32');

fwrite(fidz, max_z*10, 'float32');
fwrite(fidz, length(g(:,3)), 'float32');
fwrite(fidz, g(:,3)*10, 'float32');
fwrite(fidz2, max_z*10, 'float32');
fwrite(fidz2, length(g(:,3)), 'float32');
fwrite(fidz2, g(:,3)*10, 'float32');

fclose(fidx);fclose(fidx2);
fclose(fidy);fclose(fidy2);
fclose(fidz);fclose(fidz2);

fidrf_amp = fopen('rf_amp1801.float', 'wb');
fidrf_pha = fopen('rf_pha1801.float', 'wb');
fidrf_amp2 = fopen('rf_amp1802.float', 'wb');
fidrf_pha2 = fopen('rf_pha1802.float', 'wb');

fwrite(fidrf_amp, max_rf, 'float32');
fwrite(fidrf_amp, VoltageFactor, 'float32');
fwrite(fidrf_amp, Slice_Seperation, 'float32');
fwrite(fidrf_amp, Asymetry, 'float32');
fwrite(fidrf_amp, RewindGzMoment, 'float32');
%fwrite(fidrf_amp, zThicknessIdeal, 'float32');
%fwrite(fidrf_amp, NsliceEX, 'float32');
fwrite(fidrf_amp, length(rf), 'float32');
fwrite(fidrf_amp, abs(rf), 'float32');

fwrite(fidrf_amp2, max_rf, 'float32')
fwrite(fidrf_amp2, VoltageFactor, 'float32');
fwrite(fidrf_amp2, Slice_Seperation, 'float32');
fwrite(fidrf_amp2, Asymetry, 'float32');
fwrite(fidrf_amp2, RewindGzMoment, 'float32');
%fwrite(fidrf_amp2, zThicknessIdeal, 'float32');
%fwrite(fidrf_amp2, NsliceEX, 'float32');;
fwrite(fidrf_amp2, length(rf), 'float32');
fwrite(fidrf_amp2, abs(rf), 'float32');

fwrite(fidrf_pha, max_rf, 'float32');
fwrite(fidrf_pha, length(rf), 'float32');
fwrite(fidrf_pha2, max_rf, 'float32');
fwrite(fidrf_pha2, length(rf), 'float32');

angle_rfs = angle(factor1*rf);
angle_rfs2 = angle(factor2*rf);
index = angle_rfs<0;
angle_rfs(index) = angle_rfs(index)+2*pi*0.999;
index2 = angle_rfs2<0;
angle_rfs2(index2) = angle_rfs2(index2)+2*pi*0.999;
fwrite(fidrf_pha, angle_rfs(:), 'float32');
fwrite(fidrf_pha2, angle_rfs2(:), 'float32');

fclose(fidrf_amp);
fclose(fidrf_pha);
fclose(fidrf_amp2);
fclose(fidrf_pha2);


