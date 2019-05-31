# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
# This routine computes the Chebyshev coefficients of a function.
# Adapted from Numerical Recipes (WH Press et al., Cambridge Univ. Press, UK,
# 1989).
# In this version the function fv is given as a vector in N equidistant 
# points on the interval [a,b]. A linearly interpolated value of fv
# will be used if the point x below is not one of the sampled points.
#
#			Input:
#  tv: vector of the independent variable on the interval [a,b], where
#	N=length(tv), a=tv(1), b=tv(N);
#  fv: vector of fct. values of length N to be expanded in a series
#	in Chebyshev polynomials;
#  n: number of the Chebyshev coeff.s in the approximation.
#
#			Output:
#  ck: vector of the Chebyshev coeff.s of length n.


function ck=cheb_linip(tv,fv,n)

   N=length(tv);
   if (N!=length(fv))
     error("In cheb_linip: independent and dependent variable are\ 
            incompatible!\n");
   endif
   a=tv(1);
   b=tv(N);
   bma=0.5*(b-a);
   bpa=0.5*(b+a);
   for k=1:n
      y=cos(pi*(k-0.5)/n);
      x=y*bma+bpa;
      jn=1;
      jx=N;
      jh=floor(0.5*(jx+jn));
      while (jx-jn>1)
        if ((x==tv(jx))||(x==tv(jn))||(x==tv(jh))) break; endif
        if (x>tv(jh))
           jn=jh;
	else
	   jx=jh;
	endif
	jh=floor(0.5*(jx+jn));
      endwhile
      f(k)=fv(jn)+(fv(jx)-fv(jn))*(x-tv(jn))/(tv(jx)-tv(jn));
   endfor
   fac=2./n;
   for j=1:n
      sum=0.;
      for k=1:n
         sum=sum+f(k)*cos((pi*(j-1))*((k-0.5)/n));
      endfor
      ck(j)=fac*sum;
   endfor
   ck(1)=0.5*ck(1);

endfunction








