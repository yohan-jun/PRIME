disp('######################################');
disp('#### Design a dual density spiral ####');
disp('####                              ####');
disp('######################################');
disp(' ');



[k,g,s,time] = vdSpiralDesign(16, 0.83,[55,55,10,10],[0,0.2,0.3,1],4,15,4e-3,'cubic');

L = length(s);

figure, subplot(2,2,1), plot(k*[1;i;0]*exp(i*2*pi*[1*0])); axis([-6.5,6.5,-6.5,6.5]), title('k-space')
	subplot(2,2,2), plot(g); axis([0,L,-4.5,4.5]); title('gradient waveforms')
	subplot(2,2,3), plot(abs(g*[1;i;0])); axis([0,L,0,4.5]); title('gradient magnitude')
	subplot(2,2,4), plot(abs(s*[1;i;0])); axis([0,L,0,15.5]); title('slew-rate magnitude')

