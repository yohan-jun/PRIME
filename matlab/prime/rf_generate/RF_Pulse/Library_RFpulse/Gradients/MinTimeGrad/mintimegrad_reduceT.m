function g = mintimegrad_reduceT(Nest,g0,gf,moment,T,gmax,smax,t0,type,mintime_fac)    

T_for_opt = T*mintime_fac;
g = mintimegradOriginal(Nest,g0,gf,moment,T_for_opt,gmax,smax,0,3);
g = g.';
if size(g,1) == 3
	gx = reshape(repmat(g(1,:),mintime_fac,1),[],1);
	gy = reshape(repmat(g(2,:),mintime_fac,1),[],1);
	gz = reshape(repmat(g(3,:),mintime_fac,1),[],1);
	g = [gx, gy, gz];
elseif size(g,1) == 2
	gx = reshape(repmat(g(1,:),mintime_fac,1),[],1);
	gy = reshape(repmat(g(2,:),mintime_fac,1),[],1);
	g = [gx, gy];
elseif size(g,1) == 1
	g = reshape(repmat(g(1,:),mintime_fac,1),[],1);
end