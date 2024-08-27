function [x,info,z,w,fsbl,bndd]=socp1(f,A,b,c,d);

% [x,info,z,w,fsbl,bndd]=socp1(f,A,b,c,d);
%
% Analytic solution of the Second Order Cone Problem
% with a single constraint.
%
%  PROBLEM:
%    minimize      f'*x
%    subject to    norm(A*x+b) <= c*x+d ,  i=1,...,L
%
%   PROBLEM DESCRIPTION:
%    f -- vector defining linear objective
%    A, b, c, d -- constraint (A is m*n; b is m*1; c is 1*n; d is scalar)
%        The matrix [A;c] should be full rank, any problem can be
%        reformulated to satisfy this condition.
%
%  OUTPUT ARGUMENTS:
%    x    -- solution
%    info -- string: 'max. iterations exceeded',
%            'absolute accuracy reached', 'relative accuracy reached',
%            'target reached', 'unachievable target'
%    z, w -- solution to the dual problem
%    fsbl -- 1 if problem is feasible, 0 otherwise
%    bndd -- 1 if solution is bounded, 0 otherwise  (-1 if problem problem
%							is infeasible)

% mlobo@isl.stanford.edu -- 1997


[m,n]=size(A);
m=m+1;

eq_eps=1e-6;	% to decide when dual point is not in hyperplane (times m?)

if rank(A'*A+c'*c)<n,
	disp('Columns of [A;c] are not independent, solution not computed.');
	fsbl=-1; bndd=-1;
	info='error';
	return;
end

At=A'*A-c'*c;
bt=A'*b-c'*d;
[U,S,V]=svd(At);

if S(n,n)>n*S(1,1)*eps,  % At IS FULL RANK

	Ati=V*diag(1./diag(S))*U';
	den=f'*Ati*f;
	tp=min(min(U'*V));	% approx. -1 for hyperbolic problems,
				%     and +1 for elliptic problems

	if abs(den)<10*eps,
		if tp<0,
			disp('Degenerate hyperbolic problem.');
			x=-[A;c]\[b;d];
			if norm([A*x+b;c*x+d])<10*m*eps,
				disp('Solution is not unique.');
				fsbl=1; bndd=1;
			else
				disp(['Bounded objective, ',...
					'unbounded solution (not computed).']);
				fsbl=1; bndd=0;
			end
		else
			disp(['Cannot find solution to elliptic problem, ',...
				'(badly confitioned?)']);
			fsbl=-1; bndd=-1;
			info='error';
			return;
		end
	else
		delta=(bt'*Ati*bt-b'*b+d^2)/den;
		x=-Ati*(bt+sqrt(delta)*f);
		if ~(delta>=0 & delta<Inf & c*x+d>=0), 
			if tp<0,
				fsbl=1; bndd=0;
			else
				fsbl=0; bndd=-1;
			end
		else
			if tp>0 & delta<eps,
				disp(['Problem is marginally feasible ',...
					'(elliptic).']);
			end
			fsbl=1; bndd=1;
		end
	end

		%keyboard;

else  % At IS SINGULAR

	delta=d-b'*A*inv(A'*A)*c';
	if delta<0,  % strictly infeasible (or <-eps)
		fsbl=0; bndd=-1;
	elseif delta==0,  % marginally feasible/infeasible (or <eps)
		disp(['Problem is marginally feasible ',...
			'(parabolic).']);
		x=-[A;c]\[b;d];
		if norm([A*x+b;c*x+d])>10*m*eps,
			fsbl=0; bndd=-1;
		else
			delta=f'*inv(A'*A)*c';
			if abs(delta)<10*n*eps,
				disp('Degenerate parabolic problem, ',...
					'solution is not unique.');
				fsbl=1; fsbl=1;
			elseif delta>0
				fsbl=1; bndd=1;
			else
				fsbl=1; bndd=0;
			end
		end
	else  % strictly feasible
		V1=V(:,1:n-1);
		v2=V(:,n);
		if abs(v2'*f)<10*n*eps,  % f "horizontal"
			fsbl=1; bndd=0;
		else
			y=-[At*V1,f]\bt;	% y(n) is 1/2*lambda^(-1)
			a=bt-y(n)*f;
			if n>1,
				x1=V1*y(1:n-1);
			else
				x1=0;
			end;
			z1=(-a'*x1-b'*b+d^2)/(a'*v2);
			x=x1+v2*z1;
			if y(n)<0,
				fsbl=1; bndd=0;
			else
				fsbl=1; bndd=1;
			end
		end
	end

		%keyboard;

end


if ~fsbl,
	info='problem is infeasible';
elseif ~bndd,
	info='problem is unbounded'; 
else
	info='solution found';
end


if fsbl & bndd,  % dual solution
	u=A*x+b;
	t=c*x+d;
	k=(-A'*u+c'*t)\f;
	z=-k*u;
	w=k*t;
	err=norm(A'*z+c'*w - f);
	if err>eq_eps,
		if norm([u;t])<10*m*eps,
			z=[];
			w=[];
			disp(['Primal solution at non-differentiable ',...
				'point, dual solution not computed.']);
		else
			disp(['Warning: complementary slackness condition ',...
				'is not satisfied, ', num2str(err)]);
		end
	end
end
