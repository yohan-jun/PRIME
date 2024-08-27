%
%	An example of minimum-time VERSE applied
%	to excitation.
%

flip = 90;		% degrees.
gamma = 4258;		% gamma/2/pi at 1.5T.
T = .000004;		% Seconds.
N = 1200;		% Number of points.

x = 2*([1:N]/N-0.5);			% x from -1 to 1.
b1 = sinc(10*x);				% Sinc-RF pulse.
b1 = b1/(sum(b1)*T*gamma)*flip/360;	% Scale pulse to Gauss.

					% Note that max(b1) exceeds .14 Gauss.
					% mintverse will fix that!

g = 0.5*ones(size(b1));			% G/cm - arbitrarily chosen.

b1 = [0 b1 0];	
g = [0 g 0];


[bv,gv] = mintverse(b1,g,T,0.15,3.9,14500,T);
				% Input and output will have same step, T.

% Now, just for fun, repeat, but constrain the energy to 75% what
%	the energy of the VERSE pulse we got without an energy constraint.

ebv = sum(bv.*bv)*T;		% Energy in VERSE pulse.
[bve,gve] = mintverse(b1,g,T,0.15,3.9,14500,T, 0.75*ebv);

% Time vectors, for plotting.
tv = [1:length(bv)]' *T;
tve = [1:length(bve)]' *T;
t = [1:length(b1)]' *T;

subplot(2,3,1);
plot(t,b1);
grid on;
title('Starting B1 pulse vs t');
xlabel('Time (s)');
ylabel('B1 (G)');

subplot(2,3,2);
plot(tv,bv);
grid on;
title('VERSE B1 pulse vs t');
xlabel('Time (s)');
ylabel('B1 (G)');

subplot(2,3,3);
plot(tve,bve);
grid on;
title('Energy-Constrained VERSE B1 pulse vs t');
xlabel('Time (s)');
ylabel('B1 (G)');

subplot(2,3,4);
plot(t,g);
grid on;
title('Starting Gradient pulse vs t');
xlabel('Time (s)');
ylabel('Gradient (G/cm)');

subplot(2,3,5);
plot(tv,gv);
grid on;
title('VERSE Gradient pulse vs t');
xlabel('Time (s)');
ylabel('Gradient (G/cm)');

subplot(2,3,6);
plot(tve,gve);
grid on;
title('Energy-Constrained VERSE Gradient');
xlabel('Time (s)');
ylabel('Gradient (G/cm)');



