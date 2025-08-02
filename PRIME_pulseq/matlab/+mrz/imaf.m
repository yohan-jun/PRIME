function img = imaf(data,shift, dim, os);
% fft(shift) imab, only 2d supported so far,
% dim sos if defined
% call function  imaf(data,shift[,dim])
% data: coils * adc * phs * slc
if nargin == 1
    shift = 1;
    dim = 1;
end
if nargin < 4
    os = 2;
end
figure,
if nargin==3 || dim
    sos=@(data,dim) squeeze(sqrt(sum(abs((data).^2),dim)));
    if shift
        img = sos(ifftshift(ifftshift(ifft(ifft(data,[],2),[],3),2),3),1);
        img = img((1/2-1/os/2)*end+1:(1/2+1/os/2)*end,:,:);
        imab(img),colormap(gray),colorbar
    else
        img = sos((ifft(ifft(data,[],2),[],3)),1);
        img = img((1/2-1/os/2)*end+1:(1/2+1/os/2)*end,:,:);
        imab(img),colormap(gray),colorbar
    end
else
    if shift
        img = ifftshift(ifftshift(ifft(ifft(data,[],2),[],3),2),3);
        imab(img),colormap(gray),colorbar
    else
        img = ifft(ifft(data,[],2),[],3);
        imab(img),colormap(gray),colorbar
    end    
end


    
