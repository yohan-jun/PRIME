function [k,g,s,time] = vdSpiralDesign(Nitlv, res, fov,radius , Gmax, Smax, T, interpType)

% [k,g,s,time] = vdSpiralDesign(Nitlv, res, fov,radius , Gmax, Smax, T, interpType)
%
% This function designs a variable density spiral
%
%	Input:
%		Nitlv	-	Number of interleves
%		res	-	resolution (in mm)
%		fov	- 	vector of fov (in cm)
%		radius	-	vector of radius corresponding to the fov
%		Gmax	-	max gradient (default 4 G/CM)
%		Smax	-	max slew (default 15)
%		T	-	sampling rate (in ms)
%		interpType- 	type of interpolation used to interpolate the fov
%				accept: linear, cubic, spline
%	
%	Output:
%		k	-	the k-space trajectory
%		g	-	the gradient waveform
%		s	-	the slew rate
%		time	-	total time
%
%	example:
%		design a dual density spiral
%		
%	[k,g,s,time] = vdSpiralDesign(16, 1,[35,35,10,10],[0,0.1,0.15,1],4,15,4e-3,'linear');
%	figure, subplot(2,1,1),plot(k*[1;i;0]*exp(i*2*pi*[1:16]/16));axis('square');
%	subplot(2,1,2), plot(g), title('Gradient Waveform');
%
%
%
%	(c) Michael Lustig 2007


kmax = 5/res;

if length(radius)<2
	error(' radius must be at least length=2');
end

dr = 1/500/max(fov/Nitlv);
r = 0:dr:kmax;   kmax = max(r);

fov =  interp1(radius*kmax,fov,r,interpType);
dtheta = 2*pi*dr.*fov/Nitlv;
theta = cumsum(dtheta);

C = r.*exp(i*theta);
C = [real(C(:)), imag(C(:)), C(:)*0];
%[C,g,s,time] = minTimeGradient_mex(C, 0, Gmax, Gmax, Smax,T);
[C,time,g,s,k] = minTimeGradient(C, 0, [], Gmax, Smax,T,0);


k = C;



