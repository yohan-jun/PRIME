disp('############################################');
disp('#### Design a line (trapezoid gradient) ####');
disp('####                                    ####');
disp('############################################');
disp(' ');


kx = linspace(-5,5, 256);
%kx = linspace(0,4/0.3, 256);
ky = eps;
C = kx + i*ky;

[k,time,g,s] = minTimeGradient(C, 0, 0, 4, 15,4e-3);
L = length(s);

figure, subplot(2,2,1), plot(k); axis([-6.5,6.5,-6.5,6.5]), title('k-space')
	subplot(2,2,2), plot([real(g(:)),imag(g(:))]); axis([0,L,-4.5,4.5]); title('gradient waveforms')
	subplot(2,2,3), plot(abs(g)); axis([0,L,0,4.5]); title('gradient magnitude')
	subplot(2,2,4), plot(abs(s)); axis([0,L,0,15.5]); title('slew-rate magnitude')

