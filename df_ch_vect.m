# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
# This routine computes the derivative of a function given in n equidistant
# sampling points x over the interval [a,b].
# First, the function values are computed at the Chebyshev coeffs.s by using
# the function cheb_linip(). Then, the Chebyshev coeffs. of the derivative 
# are computed from those of the original function by means of cheb_der().
# Finally, the derivative function is evaluated at the same equidistant 
# sampling points as the original function by applying chebev_vect().
#
#			Input:
#   x: row vector of the equidistant sampling points in the interval [a,b] 
#      with x(1)=a, x(n)=b, where n=length(x);
#   f: row vector of the function values at x;
#   m: the no. of Chebyshev coeffs. of the function f.
#
#			Output:
#   f_der: row vector of the values of the derivative of f at x.
#
#			Internal variables:
#   c_f: vector of the Chebyshev coeff.s of f;
#   ch_der: vector of the Chebysev coeff.s of the derivative of f.
#
#		External fct.s and procedures:
#   cheb_linip(): computing the Chebyshev coeff.s;
#   cheb_der(): computing the Chebyshev coeff.s of the derivative of f;
#   chebev_vect(): evaluating the derivative of f at x.


function f_der=df_ch_vect(x,f,m)

	n=length(x);
	a=x(1);
	b=x(n);

#   Here the Chebyshev coeffs. of f will be computed:
	c_f=cheb_linip(x,f,m);

#   Now, the Chebyshev coeff.s of the derivative are computed from c_f:
	ch_der=cheb_der(a,b,c_f,m);

#  The derivative is evaluated at x:
	f_der=chebev_vect(a,b,ch_der,x);	# always a row vector

endfunction
