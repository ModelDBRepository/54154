# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This procedure creates the `derivative' coeff. matrix of a fct.
#	expanded in terms of Cheb. polynomials. (Cheb. approx.)
#
#			Input:
#  a,b: endpoints of the interval [a,b] where the Cheb. approx. is defined;
#  n: no. of terms in the Cheb. expansion (approx)., i.e. n-1 is the order
#	of the approx. (=highest degree among the Cheb. pol.s);
#
#			Output:
#  D: n x n matrix of the coeffs. of the Chebyshev expansions of the
#	derivatives of the Cheb. polynomials;
#  Note that the 1st column and the last row of D consist of zeros because 
#  the coeff. of T0(x) no longer appears in the derivative fct., and the
#  coeff. of Tn(x) is also zero.
#
#			External functions:
#  cheb_der(): computes the Chebyshev coeff.s of a fct. from its Chebyshev 	
		coeffs. of order n-1 on the interval [a,b].


function D=d_cch(a,b,n);

	D=zeros(n,n);
	TD=eye(n);

	for k=1:n
	   D(:,k)=cheb_der(a,b,TD(k,:),n);
	endfor

endfunction
