%WARNING: need to make sure that gradient and rf go to zero at begin and
%end of the pulse, e.g. for spiral need to have rewinder

function create_header_parallel7T(g,b_est,oversample)

B_length = length(g)*oversample;
n_coil = length(b_est)/B_length;

for count = 1:n_coil
    rf(count,:) = b_est(1+(count-1)*B_length:count*B_length);
end


%DATA READING OUT TO FILES

max_x = max(g(:,1));
max_y = max(g(:,2));
max_z = max(g(:,3));
max_rf = max(abs(rf(:)));
RFpeakVolt = max_rf

min_x = min(g(:,1));
min_y = min(g(:,2));
min_z = min(g(:,3));
%min_rf = min(rf);

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
% if abs(min_rf)> abs(max_rf)
%     max_rf = abs(min_rf);    
% end
if max_rf == 0
    max_rf = 0.1;
end

angle_rfs = angle(rf);
index = angle_rfs<0;
angle_rfs(index) = angle_rfs(index)+2*pi;
% angle_rfs(8,:).'
% cumsum(g(:,3)*10)

for (ii = 1:8)
    
    % enter tx dir
    
    dirName = ['TX' num2str(ii)];
    eval(['cd ' dirName])
  
    % do grads
    fidx = fopen('gx.float', 'wb');
    fidy = fopen('gy.float', 'wb');
    fidz = fopen('gz.float', 'wb');

    fwrite(fidx, max_x*10, 'float32');
    fwrite(fidx, length(g(:,1)), 'float32');
    fwrite(fidx, g(:,1)*10, 'float32');

    fwrite(fidy, max_y*10, 'float32');
    fwrite(fidy, length(g(:,2)), 'float32');
    fwrite(fidy, g(:,2)*10, 'float32');

    fwrite(fidz, max_z*10, 'float32');
    fwrite(fidz, length(g(:,3)), 'float32');
    fwrite(fidz, g(:,3)*10, 'float32');

    fclose(fidx);
    fclose(fidy);
    fclose(fidz);

    % do rf amp phas
    
    fname_rfamp = sprintf('rf_amp.float');
    fidrf_amp = fopen(fname_rfamp, 'wb');
    fname_rfpha = sprintf('rf_pha.float');
    fidrf_pha = fopen(fname_rfpha, 'wb');
    fwrite(fidrf_amp, max_rf, 'float32');
    fwrite(fidrf_amp, length(rf), 'float32');
    fwrite(fidrf_amp, abs(rf(ii,1:end)), 'float32');
    fwrite(fidrf_pha, max_rf, 'float32');
    fwrite(fidrf_pha, length(rf), 'float32');
    fwrite(fidrf_pha, angle_rfs(ii,1:end), 'float32');
    fclose(fidrf_amp);
    fclose(fidrf_pha);
    
    cd ..
    
end
