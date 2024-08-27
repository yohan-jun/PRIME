disp('############################################');
disp('#### Design a rosette trajectory        ####');
disp('####                                    ####');
disp('############################################');
disp(' ');

Gmx = 4;
Smx = 15;
T = 17/Gmx;
Kmx = 6;
w1 = 0.147*2*pi*Gmx;
w2 = 0.087/1.02*2*pi*Gmx;
t = 0e-3:4e-3:T;
C = Kmx*sin(w1*t).*exp(i*w2*t);



[k,time,g,s] = minTimeGradient(C, 0, 0, 4, 15,4e-3);
L = length(s);

figure, subplot(2,2,1), plot(k); axis([-6.5,6.5,-6.5,6.5]), title('k-space')
	subplot(2,2,2), plot([real(g(:)),imag(g(:))]); axis([0,L,-4.5,4.5]); title('gradient waveforms')
	subplot(2,2,3), plot(abs(g)); axis([0,L,0,4.5]); title('gradient magnitude')
	subplot(2,2,4), plot(abs(s)); axis([0,L,0,15.5]); title('slew-rate magnitude')

