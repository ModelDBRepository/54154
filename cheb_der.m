# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
# This function computes the Chebyshev coefficients of the derivative of a 
# function from given Chebyshev coefficients of the function over the interval
# [a,b].
# Adapted from Numerical Recipes (WH Press et al., Cambridge Univ. Press, UK,
# 1989).
#
#	Input:
#   a,b: endpoints of the interval [a,b];
#   ch: Chebyshev coeff.s of the (original) function;
#   n: dimension of the vector ch.
#
#	Output:
#   ch_der: vector of the Chebyshev coeff.s of the derivative.



function ch_der=cheb_der(a,b,ch,n)

	ch_der(n)=0;
	ch_der(n-1)=2*(n-1)*ch(n);
	if (n>=3)
	   n2=n-2;
	   for jj=n2:-1:1
	      ch_der(jj)=ch_der(jj+2)+2*jj*ch(jj+1);
	   endfor
	endif
	ch_der(1)=0.5*ch_der(1);
	fac=2./(b-a);
	ch_der=fac.*ch_der;

endfunction
