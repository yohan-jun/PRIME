function A = im2row(k, winSize)
%res = im2row(im, winSize)
% input k (nx ny nc) : image or kspace
% input winSize (bx by) : the size of the block
% output A ((nx-bx+1)*(ny-by+1),bx*by*nc)
% Yuxin Hu, April 22,2015
[nx,ny,nc] = size(k);
bx = winSize(1);
by = winSize(2);
A = zeros((nx-bx+1)*(ny-by+1),bx*by,nc);

count = 1;
    for y = 1 : by 
        for x = 1 : bx
            A(:,count,:) = reshape(k(x:nx-bx+x,y:ny-by+y,:),(nx-bx+1)*(ny-by+1),1,nc);
            % maybe we can first reshape "k" and then use circshift to
            % accelerate it.
            count = count + 1;
        end
    end
A = reshape(A,size(A,1),size(A,2)*size(A,3));


end

