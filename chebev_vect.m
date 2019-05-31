# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#    Vectorized evaluation of a function given as Chebyshev approx.:
#		f(x)=sum[c_j*T_(j-1)(x)]
#    where c_j are the Chebyshev coefficients, and x is a vector.
#    This subroutine creates matrix A that contains the values of
#    the Chebyshev polynomials at the different components of x, i.e.
#		A(j,i)=T_(j-1)(x_i),  j=1,...,m
#    Then the function value, as vector corresponding to x is given as
#		fx=c*A
#    where c is the vector of the Chebyshev coefficients of T_(j-1)(x).
#
#		Input:
#  a: lower endpoint of the interval in which x is given;
#  b: upper endpoint of the interval in which x is given;
#  c: vector of the Chebyshev coefficients (see above);
#  x: vector at whose components the Chebyshev approx. of f(x) is 
#	to be evaluated.
#
#		Output:
#   fx: row vector of the function values at x.



function fx=chebev_vect(a,b,c,x)

	m=length(c);
	if (m<1)
	   error("No Chebyshev coeffs.!\n");
	endif
	if (rows(c)>1) c=c'; endif

	n=length(x);
	if (rows(x)>1) x=x'; endif


	y=(2*x-a-b)/(b-a);	# interval transformation
	y2=2*y;

	A=zeros(m,n);
	A(1,:)=ones(1,n);

	if(m==1) 
	   fx=c*A;
	   return   
	endif

	A(2,:)=y;
	if (m==2) 
	   fx=c*A;
	   return  
	endif

	for jj=3:m
	   A(jj,:)=y2.*A(jj-1,:)-A(jj-2,:);
	endfor

	fx=c*A;

endfunction
