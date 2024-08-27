% design a waveform for spiral excitation.
%[k,g,s,time] = vdSpiralDesign(1, 10, [25,25],[0,1] , 4, 15, 4e-3, 'linear');
%g = g(end:-1:1,[1,2]);
%s = s(end:-1:1,[1,2]);
%figure, subplot(2,2,1), plot(k*[1;i;0]);
%	subplot(2,2,2), plot(g);
%	subplot(2,2,3), plot(abs(g*[1;i]));
%	subplot(2,2,4), plot(abs(s*[1;i]));


% design multiple density spiral
[k,g,s,time] = vdSpiralDesign(16, 0.83,[55,55,10,10],[0,0.2,0.3,1],4,15,4e-3,'cubic');
g = g(:,[1,2]);
s = s(:,[1,2]);
L = length(s);
figure, subplot(2,2,1), plot(k*[1;i;0]*exp(i*2*pi*[1*0])); axis([-6.5,6.5,-6.5,6.5])
	subplot(2,2,2), plot(g); axis([0,L,-4.5,4.5]);
	subplot(2,2,3), plot(abs(g*[1;i])); axis([0,L,0,4.5]);
	subplot(2,2,4), plot(abs(s*[1;i])); axis([0,L,0,15.5]);


% design a circular trajectory
C = exp(i*2*pi*linspace(0,1,512))*6;
[k,time,g,s] = minTimeGradient(C, [], [], 4, 15,4e-3);
L = length(s);
figure, subplot(2,2,1), plot(k);axis([-6.5,6.5,-6.5,6.5])
	subplot(2,2,2), plot([real(g(:)),imag(g(:))]);  axis([0,L,-4.5,4.5]);
	subplot(2,2,3), plot(abs(g));axis([0,L,0,4.5]);
	subplot(2,2,4), plot(abs(s)); axis([0,L,0,15.5]);



% design a line

%kx = linspace(-5,5,256); kx = [kx, -kx, kx, -kx, kx, -kx];
kx = linspace(-5,5, 256);
ky = eps;
C = kx + i*ky;
%C = filter2(ones(1,15)/15,C,'crop');
[k,time,g,s] = minTimeGradient(C, 0, 0, 4, 15,4e-3);
L = length(s);
figure, subplot(2,2,1), plot(k);axis([-5.5,5.5,-5.5,5.5])
	subplot(2,2,2), plot([real(g(:)),imag(g(:))]); axis([0,L,-4.5,4.5]);
	subplot(2,2,3), plot(abs(g));axis([0,L,0,4.5]);
	subplot(2,2,4), plot(abs(s)); axis([0,L,0,15.5]);




% design a Rosette trajectory
Gmx = 4;
Smx = 15;
T = 17/Gmx;
Kmx = 6;
%w2 = 2*pi/17e-3;
w1 = 0.147*2*pi*Gmx;
w2 = 0.087/1.02*2*pi*Gmx;
t = 0e-3:4e-3:T;
G = 1/4.257*Kmx/2*((w1+w2)*exp(i*(w1+w2)*t) + (w1-w2)*exp(-i*(w1-w2)*t));
S = 1/4.257*Kmx/2*i*((w1+w2)^2*exp(i*(w1+w2)*t) -(w1-w2)^2*exp(-i*(w1-w2)*t));

dt = [100, exp(-linspace(0,1,length(G)-1)*300)*20+1];
tt = cumtrapz(dt*4e-3);
newt = 0:4e-3:tt(end);
K = cumtrapz(G*4e-3*4.257);
KK = interp1(tt,K,newt,'spline');
GG = diff(KK/4e-3/4.257);


C = Kmx*sin(w1*t).*exp(i*w2*t);


[k,time,g,s] = minTimeGradient(C, 0, [], Gmx, Smx,4e-3);
L = length(s);
figure, subplot(2,2,1), plot(k);axis([-6.5,6.5,-6.5,6.5])
	subplot(2,2,2), plot([real(g(:)),imag(g(:))]);  axis([0,L,-4.5,4.5]);
	subplot(2,2,3), plot(abs(g));axis([0,L,0,4.5]);
	subplot(2,2,4), plot(abs(s));axis([0,L,0,15.5]);




