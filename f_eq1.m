# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
# This is a nonlinear equation used in the estimation of the par. values
# in alpha(V)=ap*x/(exp(x)-1), where x=q*V+b.

function y=f_eq1(x)

	global ap ri;

	y(1)=exp(x)-1-ap*x/ri;

endfunction
