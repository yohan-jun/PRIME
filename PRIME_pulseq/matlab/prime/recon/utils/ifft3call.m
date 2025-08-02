function res = ifft3call(x)
fctr = size(x,1)*size(x,2)*size(x,3);

X = fftshift(ifft(ifftshift(x,1),[],1),1);

X = fftshift(ifft(ifftshift(X,2),[],2),2);

res = fftshift(ifft(ifftshift(X,3),[],3),3) * sqrt(fctr);



