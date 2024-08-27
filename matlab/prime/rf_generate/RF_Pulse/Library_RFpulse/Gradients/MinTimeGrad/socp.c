/*
 *   socp.c
 *
 *   Second-Order Cone Programming
 *   C implementation
 *
 *   mlobo@isl.stanford.edu -- 96/97
 */

#include <math.h>
#include <float.h>

#include "socp.h"



void dzero(
 int n,
 double *p
)
     /*
      *  set vector of n doubles to zero
      */
{
  /* faster version:
    n*=sizeof(double)/sizeof(char);
    memset(p,0,n);
  */
  /* alternative: (using BLAS to avoid need for memory.h) */
  int int1=1;
  double double0=0.0;
  dscal_(&n,&double0,p,&int1);
}



double dsumdiv(
 int n,
 double *x,
 double *y
)
     /*
      * returns sum of elementwise division of x by y
      *   in matlab notation: sum(x./y)
      */
{
  int i;
  double sumdiv=0.0;
  for (i=0; i<n; ++i)
    sumdiv += x[i]/y[i];
  return sumdiv;
}



void dupge(
 int n,
 double *A
)
     /*
      *  Convert symmetric matrix from upper storage to full storage
      */
{
  int i, j, k;
  int int1=1;

  for (i=n-1, j=1, k=n; i>0; --i, j+=n+1, k+=n+1)
    dcopy_(&i,A+k,&n,A+j,&int1);
}



void dgapdev(
 int m,
 int L,
 int *N,
 double *u,
 double *z,
 double *pgap,
 double *pdev
)
     /*
      *  note: this fct currently not used
      *
      *  compute duality gap and deviation from centrality
      *
      *  output:
      *   *pgap = duality gap
      *   *pdev = deviation from centrality
      */
{
  int i, j, k;
  double fu, fz;

  int int1=1;

  /* gap = u'*z; */
  *pgap=ddot_(&m,u,&int1,z,&int1);

  /* dev = -sum([log(SF'*(u.^2));log(SF'*(z.^2))]) + 2*L*(log(gap)-log(L)); */
  *pdev=2*L*(log(*pgap)-log(L));
  for (i=0, k=0; i<L; ++i) {
    for (j=0, fu=0.0, fz=0.0; j<N[i]-1; ++j, ++k) {
      fu-=SQR(u[k]);
      fz-=SQR(z[k]);
    }
    fu+=SQR(u[k]);
    fz+=SQR(z[k]); ++k;
    *pdev-=log(fu)+log(fz);
  }
}



void dgrad(
   /* input args. */
 double w,
 int L,
 double c1,      /* pre-computed constants, dependent on: u, z, du, dz */
 double c2,
 double c3,
 double *d1,
 double *d2,
 double *d3,
 double *e1,
 double *e2,
 double *e3,
 double p,       /* du scaling */
 double q,       /* dz scaling */
   /* output args. */
 double *pgp,    /* return gradient */
 double *pgq,
 double *t1,     /* intermediate values that will be re-used in dhess() */
 double *t2,
 double *t3,
 double *t4
)
     /*
      *  compute gradient of potential fct wrt. p and q
      *   at u+p*du, z+q*dz
      *
      *  output:
      *   *pgp = gradient wrt. p
      *   *pgq = gradient wrt. q
      */
{
  double ptwo=2*p;
  double psqr=SQR(p);
  double qtwo=2*q;
  double qsqr=SQR(q);

  int int1=1;

  /* t1 = d2 + p*d3 */
  dcopy_(&L,d2,&int1,t1,&int1);
  daxpy_(&L,&p,d3,&int1,t1,&int1);

  /* t2 = d1 + 2*p*d2 + p*p*d3 */
  dcopy_(&L,d1,&int1,t2,&int1);
  daxpy_(&L,&ptwo,d2,&int1,t2,&int1);
  daxpy_(&L,&psqr,d3,&int1,t2,&int1);

  /* t3 = e2 + q*e3 */
  dcopy_(&L,e2,&int1,t3,&int1);
  daxpy_(&L,&q,e3,&int1,t3,&int1);

  /* t4 = e1 + 2*q*e2 + q*q*e3 */
  dcopy_(&L,e1,&int1,t4,&int1);
  daxpy_(&L,&qtwo,e2,&int1,t4,&int1);
  daxpy_(&L,&qsqr,e3,&int1,t4,&int1);

  *pgp=w*c2/(c1+p*c2+q*c3)-2*dsumdiv(L,t1,t2);
  *pgq=w*c3/(c1+p*c2+q*c3)-2*dsumdiv(L,t3,t4);
}



void dhess(
   /* input args. */
 int L,
 double *d3,      /* pre-computed constants */
 double *e3,
 double *t1,      /* intermediate values from dgrad() */
 double *t2,
 double *t3,
 double *t4,
   /* output args. */
 double *php,
 double *phq
)
     /*
      *  compute hessian of primal barrier wrt. p and q
      *   at u+p*du, z+q*dz
      *
      *  output:
      *   *php = hessian wrt. p
      *   *phq = hessian wrt. q
      */
{
  int i;

  for (i=0, *php=0.0, *phq=0.0; i<L; ++i) {
    *php+=(d3[i]*t2[i]-2*SQR(t1[i]))/SQR(t2[i]);
    *phq+=(e3[i]*t4[i]-2*SQR(t3[i]))/SQR(t4[i]);
  }
  *php*=-2;
  *phq*=-2;
}



void socp_getwork(
	/* input args.: */
 int L,
 int *N,
 int n,
 int max_iter,
 int out_mode,
        /* output args.: */
 int *mhist,
 int *nhist,
 int *ndbl,
 int *nint
)
     /*
      *  returns workspace size in number of doubles, pointers and ints
      *  required by socp()
      *  (use a macro in socp.h instead?)
      *
      *  input arguments:
      *   use the same values as for socp()
      *
      *  output arguments:
      *   *mhist-- number of lines in *hist (elements are doubles)
      *   *nhist-- number of columns in *hist
      *   *ndbl -- number of doubles in *dblwork
      *   *nptr -- number of pointers in *ptrwork
      *   *nint -- number of integers in *intwork
      */
{
  int i;
  int m;

  /* compute m */
  for (i=0, m=0; i<L; m+=N[i++]);

  *mhist = out_mode;   /* 0: none;  1: gap only;  2: gap and deviation */
  *nhist = max_iter+1; /* for initial point and after each iteration */

  *ndbl = 7*m + 2*n + 10*n + 2*SQR(n) + 11*L;
  *nint = n;
}



int socp(
 int L,
 int *N,
 int n,

 double *f,
 double *A,
 double *b,

 double *x,
 double *z,

 double abs_tol,
 double rel_tol,
 double target,
 int *iter,

 double Nu,

 int *info,
 int out_mode,
 double *hist,

 double *dblwork,
 int *intwork
)
     /*
      *  solve second order cone problem
      *
      */
{
  int m;
  double w;
  double gap;
  double dev;

  int iter_out;
  double lbound;
  double *u, *fc, *gup, *gu;    /* u[m], fc[L], gup[m], gu[m] */
  double *du, *dz;              /* du[m], dz[m] */
  double *gx, *dx;	        /* gx[n], dx[n] */
  double *Hx;                   /* Hx[n*(n+1)] */
  double s;

    /* to record reason for exiting main loop */
  int xabs=0;
  int xrel=0;
  int xtarget=0;
  int xiter=0;

    /* for plane search */
  double c1, c2, c3;
  double *d1, *d2, *d3, *e1, *e2, *e3;    /* [L] */
  double p, q, gp, gq, hp, hq, dp, dq;
  double *t1, *t2, *t3, *t4;              /* [L] */
  double lambda2;
  int iter_plane;

    /* for line search */
  double p0, q0;
  double alpha, galpha;
  double *ua, *za;              /* ua[m], za[m] */
  double fu, fz;
  int out_barrier;

    /* general */
  int i, j, k, l;

    /* for linear system solver */
    /* for both dposvx_() and dgelss_() */
  int use_ls = 0;     /* will be changed to 1 if dgelss_ is to be used */
  int la_info;        /* exit info from dposvx_ and dgelss_ */
  double *Hwork;      /* Hwork[SQR(n)] (dgelss_ only needs n) */
  double rcond;
    /* for dposvx_() */
  int equed;
  double scale, ferr, berr;
    /* for dgelss_() */
  int rank;
  int lwork = 10*n;   /* size of workspace, minimum is 5*n */

    /* constant variables for BLAS */
  int int1=1;
  double double1=1.0, double0=0.0, double_1=-1.0;



  /* compute m: size of dual variable */
  for (i=0, m=0; i<L; m+=N[i++]);

  /*
   * organize workspace
   *
   * total space needed =
   *   doubles:   7*m + 2*n + max(3*n,lwork) + 2*SQR(n) + 11*L
   *   ints:      n
   */

      /* dblwork is used for dposvx_() and dgelss_(), */
      /*   need 3*n and lwork doubles respectively. */
  u   = dblwork + lwork;    /* lwork=10*n */
  fc  = u + m;
  gup = fc + L;
  gu  = gup + m;
  du  = gu + m;
  dz  = du + m;
  gx  = dz + m,
  dx  = gx + n;
  Hx  = dx + n;
  d1  = Hx + SQR(n);
  d2  = d1 + L;
  d3  = d2 + L;
  e1  = d3 + L;
  e2  = e1 + L;
  e3  = e2 + L;
  t1  = e3 + L;
  t2  = t1 + L;
  t3  = t2 + L;
  t4  = t3 + L;
  ua  = t4 + L;
  za  = ua + m;
  Hwork = za + m;
  /* Hwork needs SQR(n) doubles */


  /* gap reduction vs. centering */
  w=2*L+Nu*sqrt(2*L);

  /* u=A*x+b */
  dcopy_(&m,b,&int1,u,&int1);  /* u=b */
  dgemv_("N",&m,&n,&double1,A,&m,x,&int1,&double1,u,&int1);  /* u=A*x+u */

  /* compute gap (and deviation from centrality and store in hist) */
  if (out_mode == 2) {
    dgapdev(m,L,N,u,z,&gap,&dev);
    hist[0]=gap;
    hist[1]=dev;
  } else {
    gap=ddot_(&m,u,&int1,z,&int1);  /* gap = u'*z; */
    if (out_mode == 1)
      hist[0]=gap;
  }

  /* outer loop */
  iter_out=0;
  while (!(( rel_tol<0.0
	     && (xtarget = (ddot_(&n,f,&int1,x,&int1)<target
			 || -ddot_(&m,b,&int1,z,&int1)>=target)))
	   || ( abs_tol>0.0
		&& (xabs = (gap<=abs_tol)))
	   || ( rel_tol>0.0
		&& (xrel = (((lbound=-ddot_(&m,b,&int1,z,&int1))>0.0
				&& gap/lbound<=rel_tol)
			       ||((lbound=-ddot_(&n,f,&int1,x,&int1))>0.0
				  && gap/lbound<=rel_tol))))
	   || (xiter = (iter_out>=*iter)))) {
    ++iter_out;

    /* compute gup (gradient of primal barrier wrt. u) */
    /* also, compute fc(i)=2/(t(i)^2-u(i)'*u(i)) */
    for (i=0, k=0; i<L; ++i) {
      for (j=0, fc[i]=0.0; j<N[i]-1; ++j)
	fc[i]-=SQR(u[k+j]);
      fc[i]+=SQR(u[k+j]);
      fc[i]=2.0/fc[i];
      for (j=0; j<N[i]-1; ++j, ++k)
	gup[k]=fc[i]*u[k];
      gup[k]=-fc[i]*u[k]; ++k;
    }

    /* compute gu (gradient of potential wrt. u) */
    s=w/gap;
    dcopy_(&m,gup,&int1,gu,&int1);   /* gu=gup */
    daxpy_(&m,&s,z,&int1,gu,&int1);  /* gu=gu+(w/gap)*z */

    /* compute Hx = A'*Hu*A */
    /* where   Hu(i) = fc(i)*diag([1 1 ... 1 -1]) + gup(i)*gup(i)' */
    /* Hx=0 */
    dzero(SQR(n),Hx);

    for (i=0, k=0; i<L; k+=N[i++]) {    /* for each constraint */
      /* n. of rows of A(i) */
      j = N[i]-1;

      /* Hx = Hx + fc(i)*A(i)'*A(i) */
      if (j>0)  dsyrk_("U","T",&n,&j,&(fc[i]),&(A[k]),&m,&double1,Hx,&n);

      /* Hx = Hx - fc(i)*c(i)*c(i)'  (rank one update) */
      s = -fc[i];
      dsyr_("U",&n,&s,&(A[k+j]),&m,Hx,&n);

      /* gx = [A(i);c(i)']'*gup(i)  (this is not gx, just used as aux.) */
      dgemv_("T",&(N[i]),&n,&double1,&(A[k]),&m,
	     &(gup[k]),&int1,&double0,gx,&int1);

      /* Hx = Hx + gx*gx'  (rank one update) */
      dsyr_("U",&n,&double1,gx,&int1,Hx,&n);
    }

    /* solve linear system: dx = -Hx\(A'*gu) */
    /* gx = -A'*gu */
    dgemv_("T",&m,&n,&double_1,A,&m,gu,&int1,&double0,gx,&int1);
    /* dx = Hx\gx */
    if (!use_ls) {   /* solve linear system by QR fact. */
      dposvx_("N","U",&n,&int1,Hx,&n,Hwork,&n,&equed,&scale,
	      gx,&n,dx,&n,
	      &rcond,&ferr,&berr,dblwork,intwork,&la_info);
      if (la_info>0) /* from 1 to n, Hessian not positive def.; */
	use_ls = 1;  /* n+1, Hessian badly conditioned; */
                     /* do SVD now and switch to SVD for all */
                     /* future iterations */
    }
    if (use_ls) {   /* solve linear system in least squares sense using SVD */
      dupge(n,Hx);  /* convert to general storage */
      rcond=-1;     /* keep singular values down to machine precision */
      dgelss_(&n,&n,&int1,Hx,&n,
	      gx,&n,
	      Hwork,&rcond,&rank,   /* (only first n of Hwork used, for S) */
	      dblwork,&lwork,&la_info); /* dblwork: lwork doubles (>= 5*n) */
      dcopy_(&n,gx,&int1,dx,&int1); /* dx=gx (dgelss_ overwrites gx) */
    }
    if (la_info)
      return la_info;    /* abort: return error in lapack solver */

    /* du = A*dx */
    dgemv_("N",&m,&n,&double1,A,&m,dx,&int1,&double0,du,&int1);

    /* dz = -(gu+Hu*du) */
    /* computed one constraint at a time: */
    /* dz(i)= -(gu(i)+fc(i)*[du(i);-dt(i)]+gup(i)*(gup(i)'*[du(i);dt(i)])) */
    /* dz = gu */
    dcopy_(&m,gu,&int1,dz,&int1);
    for (i=0, k=0; i<L; k+=N[i++]) {
      j = N[i]-1;
      /* dz(i) = dz(i) + fc[i]*[du(i);-dt(i)] */
      daxpy_(&j,&(fc[i]),&(du[k]),&int1,&(dz[k]),&int1);
      dz[k+j] -= fc[i]*du[k+j];
      /* s = gup(i)'*du(i) */
      s = ddot_(&(N[i]),&(gup[k]),&int1,&(du[k]),&int1);
      /* dz(i) = dz(i) + s*gup(i) */
      daxpy_(&(N[i]),&s,&(gup[k]),&int1,&(dz[k]),&int1);
    }
    /* dz = - dz */
    dscal_(&m,&double_1,dz,&int1);
    /* optional: scale dz by 1/rho=gap/w
      s=gap/w;
      dscal_(&m,&s,dz,&int1);
    */

    /*
     *  constants for plane search
     */
    c1=gap;
    c2=ddot_(&m,du,&int1,z,&int1);
    c3=ddot_(&m,u,&int1,dz,&int1);
    for (i=0, k=0; i<L; ++i) {
      d1[i]=0.0;
      d2[i]=0.0;
      d3[i]=0.0;
      e1[i]=0.0;
      e2[i]=0.0;
      e3[i]=0.0;
      for (j=0; j<N[i]-1; ++j, ++k) {
	d1[i]-=SQR(u[k]);
	d2[i]-=u[k]*du[k];
	d3[i]-=SQR(du[k]);
	e1[i]-=SQR(z[k]);
	e2[i]-=z[k]*dz[k];
	e3[i]-=SQR(dz[k]);
      }
      d1[i]+=SQR(u[k]);
      d2[i]+=u[k]*du[k];
      d3[i]+=SQR(du[k]);
      e1[i]+=SQR(z[k]);
      e2[i]+=z[k]*dz[k];
      e3[i]+=SQR(dz[k]); ++k;
    }

    /* plane search loop */
    p=0.0;
    q=0.0;
    for (lambda2=2*MAX_LAMBDA2, iter_plane=0;
	 lambda2>MAX_LAMBDA2 && iter_plane<MAX_ITER_PLANE;
	 ++iter_plane) {

      /* compute gradient and Hessian wrt. p and q */
      /*   at u+p*du, z+q*dz */
      dgrad(w,L,c1,c2,c3,d1,d2,d3,e1,e2,e3,p,q,&gp,&gq,t1,t2,t3,t4);
      dhess(L,d3,e3,t1,t2,t3,t4,&hp,&hq);

      /* Newton step */
      dp=-gp/hp;
      dq=-gq/hq;

      /* line search loop: scale down step */
      /*   until inside feasible region */
      /*   (and until between previous point and line minimum) */
      alpha=1;  /* scaling factor for dp, dq */
      p0=p;
      q0=q;
      while (1) {
	p=p0+alpha*dp;
	q=q0+alpha*dq;

	/* ua=u+p*du */
	dcopy_(&m,u,&int1,ua,&int1);
	daxpy_(&m,&p,du,&int1,ua,&int1);
	/* za=z+q*dz */
	dcopy_(&m,z,&int1,za,&int1);
	daxpy_(&m,&q,dz,&int1,za,&int1);

	/* check constraints: */
	/*   fu = t(i)^2 - u(i)'*u(i) > 0  and  t(i) > 0 */
	/*   fz = s(i)^2 - z(i)'*z(i) > 0  and  s(i) > 0 */
	for (i=0, k=0, out_barrier=0; i<L && !out_barrier; ++i) {
	  for (j=0, fu=0.0, fz=0.0; j<N[i]-1; ++j, ++k) {
	    fu-=SQR(ua[k]);
	    fz-=SQR(za[k]);
	  }
	  fu+=SQR(ua[k]);
	  fz+=SQR(za[k]);
	  if (fu<=0 || fz<=0 || ua[k]<=0 || za[k]<=0) {
	    out_barrier=1;                    /* (replace <=0 with <TINY ?) */
	    break;
	  }
	  ++k;
	}

	if (!out_barrier) {
	  /* compute gradient along search line wrt. alpha */
	  dgrad(w,L,c1,c2,c3,d1,d2,d3,e1,e2,e3,p,q,&gp,&gq,t1,t2,t3,t4);
	  galpha=dp*gp+dq*gq;
	  /* exit if: all barriers ok (i.e. out_barrier==0) */
	  /*  and gradient negative (=> between initial and minimum) */
	  if (galpha<=0)
	    break;    /* EXIT line search */
	}
	alpha/=DIV_ALPHA;         /* scale down the p, q step */
	if (alpha<MIN_ALPHA) {    /* line search failed */
	  alpha=0.0;
	  p=p0;
	  q=q0;
	  break;    /* EXIT line search */
	}
      }    /* end of line search loop */

      /* plane search convergence criterium */
      lambda2=hp*SQR(dp)+hq*SQR(dq);

    }    /* end of plane search loop */

    /* x=x+p*dx */
    daxpy_(&n,&p,dx,&int1,x,&int1);
    /* z=z+q*dz */
    daxpy_(&m,&q,dz,&int1,z,&int1);
    /* u=A*x+b*/
    dcopy_(&m,b,&int1,u,&int1);
    dgemv_("N",&m,&n,&double1,A,&m,x,&int1,&double1,u,&int1);

    /* update gap (and deviation from centrality and store in hist) */
    if (out_mode == 2) {
      dgapdev(m,L,N,u,z,&gap,&dev);
      hist[2*iter_out]=gap;
      hist[2*iter_out+1]=dev;
    } else {
      gap=ddot_(&m,u,&int1,z,&int1);  /* gap = u'*z; */
      if (out_mode == 1)
	hist[iter_out]=gap;
    }
  }    /* end of outer loop */
  *iter=iter_out;

    /* report reason for exit */
  *info = xabs + 2*xrel + 3*xtarget + 4*xiter;
    /* normal exit */
  return 0;
}    /* end of socp() */
