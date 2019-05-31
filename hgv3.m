# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This procedure computes the matrix H of g*V=Hv
#
#		Input:
#  N: order of the Cheb. approx. of the (unknown) fct. V 
#	with Cheb. coeffs v(k), k=1,...,N+1;
#  g: vector of Cheb.coeffs.of a (known) fct. dim g=n;

#		Output:
#  H: (N+n)x(N+1) matrix of the coeffs. of v(k)


function H=hgv3(N,g)
   n=length(g);
   N2=N+n;
   n=n-1;
   H=zeros(N2,N+1);
   for k=0:N
      k1=k+1;
      for l=0:n
	ah=g(l+1);
	s1=k+l+1;
	s2=k-l;
	if (s2<0) s2=-s2; endif
	s2=s2+1;
	H(s1,k1)=H(s1,k1)+ah;
	H(s2,k1)=H(s2,k1)+ah;
      endfor
   endfor
   H=0.5*H;
endfunction
