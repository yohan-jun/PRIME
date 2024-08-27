function [angle1, angle2] = abr3D(Bx,By,gx,gy,gz,x,y,z,t)

% Bx is the excitation field (ie rf) : is a vector, By is used if the rf is
% not in phase, o.w. it is set to a zero vector
% gx is how the x-direction gradient changes with time :is a vector , similar for gy and gz
% x,y,z are the vector specifying the range 
% t is the time length of the pulse in ms

l = length(Bx);
timegap = t/l;
gamma = 2*pi*4.258;
constant1 = gamma*timegap;

%initialisation
A = ones(length(y),length(x),length(z));
I = A;
B = zeros(length(y),length(x),length(z));
X_cube = B;Y_cube = B;Z_cube = B; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for count = 1:length(x)
    X_cube(:,count,:) = x(count);
end
for count = 1:length(y)
    Y_cube(count,:,:) = y(count);
end
for count = 1:length(z)
    Z_cube(:,:,count) = z(count);
end

for count = 1:l
    G = gx(count)*X_cube + gy(count)*Y_cube + gz(count)*Z_cube;
    norm = sqrt( Bx(count)^2*I + By(count)^2*I+ G.^2);
    angle = -constant1*norm;
    index1 = angle==0; norm(index1) =1;%dummie value so will not be dividing by zero when calculate direction
    index2 = angle~=0;
    direction1 = Bx(count)*I./norm;
    direction2 = By(count)*I./norm;
    direction3 = G./norm;
    alpha = cos(angle/2) - i*direction3.*sin(angle/2);
    beta = -i*(direction1 + i*direction2).*sin(angle/2);
    Anew = alpha.*A - conj(beta).*B;
    Bnew = beta.*A + conj(alpha).*B;
    A(index2) = Anew(index2); B(index2) = Bnew(index2);                
end
angle1 = A;
angle2 = B;
