function result = LORAKS_operators(x, N1, N2, Nc, R, LORAKS_type, weights)
if LORAKS_type == 1     % S matrix
    x = reshape(x,N1*N2,Nc);

    [in1,in2] = meshgrid(-R:R,-R:R);
    i = find(in1.^2+in2.^2<=R^2);

    in1 = in1(i)';
    in2 = in2(i)';

    patchSize = numel(in1);

    result = zeros(patchSize,2, Nc, (N1-2*R-even(N1))*(N2-2*R-even(N2))*2,'like',x);

    k = 0;
    for i = R+1+even(N1):N1-R
        for j = R+1+even(N2):N2-R
            k = k+1;
            Ind = sub2ind([N1,N2],i+in1,j+in2);
            Indp = sub2ind([N1,N2],-i+in1+2*ceil((N1-1)/2)+2,-j+in2+2*ceil((N2-1)/2)+2);

            tmp = x(Ind,:)-x(Indp,:);
            result(:,1,:,k) = real(tmp);
            result(:,2,:,k) = -imag(tmp);

            tmp = x(Ind,:)+x(Indp,:);
            result(:,1,:,k+end/2) = imag(tmp);
            result(:,2,:,k+end/2) = real(tmp);
        end
    end

    result = reshape(result, patchSize*Nc*2,(N1-2*R-even(N1))*(N2-2*R-even(N2))*2);

elseif LORAKS_type == -1 % S^H (Hermitian adjoint of of S matrix) 
    result = zeros(N1*N2,Nc,'like',x);

    [in1,in2] = meshgrid(-R:R,-R:R);
    i = find(in1.^2+in2.^2<=R^2);

    in1 = in1(i)';
    in2 = in2(i)';
    patchSize = numel(in1);
    nPatch = (N1-2*R-even(N1))*(N2-2*R-even(N2));

    x = reshape(x, [patchSize*2, Nc, (N1-2*R-even(N1))*(N2-2*R-even(N2))*2]);

    k = 0;
    for i = R+1+even(N1):N1-R
        for j = R+1+even(N2):N2-R
            k = k+1;
            Ind = sub2ind([N1,N2],i+in1,j+in2);
            Indp = sub2ind([N1,N2],-i+in1+2*ceil((N1-1)/2)+2,-j+in2+2*ceil((N2-1)/2)+2);

            result(Ind,:) = result(Ind,:) + complex(x(1:patchSize,:,k) + x(patchSize+1:2*patchSize,:,nPatch+k), x(1:patchSize,:,nPatch+k) - x(patchSize+1:2*patchSize,:,k)); 
            result(Indp,:) = result(Indp,:) + complex( - x(1:patchSize,:,k) + x(patchSize+1:2*patchSize,:,nPatch+k), x(1:patchSize,:,nPatch+k) + x(patchSize+1:2*patchSize,:,k)); 

        end
    end

    result = vect(result);

elseif LORAKS_type == 2  % C matrix
    x = reshape(x,N1*N2,Nc);

    [in1,in2] = meshgrid(-R:R,-R:R);
    i = find(in1.^2+in2.^2<=R^2);

    in1 = in1(i)';
    in2 = in2(i)';

    patchSize = numel(in1);

    result = zeros(patchSize*Nc,(N1-2*R-even(N1))*(N2-2*R-even(N2)),'like',x);

    k = 0;
    for i = R+1+even(N1):N1-R
        for j = R+1+even(N2):N2-R
            k = k+1;
            ind = sub2ind([N1,N2],i+in1,j+in2);
            result(:,k) = vect(x(ind,:));
        end
    end

elseif LORAKS_type == -2 % C^H (Hermitian adjoint of of C matrix) 
    result = zeros(N1*N2,Nc,'like',x);

    [in1,in2] = meshgrid(-R:R,-R:R);
    i = find(in1.^2+in2.^2<=R^2);

    in1 = in1(i)';
    in2 = in2(i)';
    patchSize = numel(in1);

    k = 0;
    for i = R+1+even(N1):N1-R
        for j = R+1+even(N2):N2-R
            k = k+1;
            ind = sub2ind([N1,N2],i+in1,j+in2);
            result(ind,:) = result(ind,:)+ reshape(x(:,k),patchSize,Nc);
        end
    end
    result = vect(result);

elseif LORAKS_type == 3  % W matrix
    W1 = weights(:,:,1:Nc,1);
    W2 = weights(:,:,1:Nc,2);
    
    W1x = reshape(W1(:).*x(:),N1*N2,Nc);
    W2x = reshape(W2(:).*x(:),N1*N2,Nc);
    
    [in1,in2] = meshgrid(-R:R,-R:R);
    i = find(in1.^2+in2.^2<=R^2);

    in1 = in1(i)';
    in2 = in2(i)';

    patchSize = numel(in1);
    nPatch = (N1-2*R-even(N1))*(N2-2*R-even(N2));
    result = zeros(patchSize*Nc,(N1-2*R-even(N1))*(N2-2*R-even(N2))*2,'like',x);

    k = 0;
    for i = R+1+even(N1):N1-R
        for j = R+1+even(N2):N2-R
            k = k+1;
            ind = sub2ind([N1,N2],i+in1,j+in2);
            result(:,k) = vect(W1x(ind,:));
            result(:,nPatch+k) = vect(W2x(ind,:));
        end
    end
    
elseif LORAKS_type == -3  % W^H (Hermitian adjoint of of W matrix)
    W1 = reshape(weights(:,:,1:Nc,1),[N1*N2,Nc]);
    W2 = reshape(weights(:,:,1:Nc,2),[N1*N2,Nc]);
    
    result1 = zeros(N1*N2,Nc,'like',x);
    result2 = zeros(N1*N2,Nc,'like',x);
    
    [in1,in2] = meshgrid(-R:R,-R:R);
    i = find(in1.^2+in2.^2<=R^2);

    in1 = in1(i)';
    in2 = in2(i)';
    patchSize = numel(in1);
    nPatch = (N1-2*R-even(N1))*(N2-2*R-even(N2));

    k = 0;
    for i = R+1+even(N1):N1-R
        for j = R+1+even(N2):N2-R
            k = k+1;
            ind = sub2ind([N1,N2],i+in1,j+in2);
            result1(ind,:) = result1(ind,:)+ reshape(x(:,k),patchSize,Nc);
            result2(ind,:) = result2(ind,:)+ reshape(x(:,nPatch+k),patchSize,Nc);
        end
    end
    result = vect(conj(W1).*result1 + conj(W2).*result2);
end
end