# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#	Steady-state current balance equation for computing Er, the
#	resting potential of a neurone:
#	   IT+IL=0
#	The function fEr2 below is used by fsolve a solver routine of
#	systems of nonlinear equations.
#	The parameters are passed as global var.s to the function.
#	The fct. fmh() computes the time course of m(t), h(t) or m(t)^k*h(t)
#	   of a membrane current when t, V(t) and init value(s) are given;



function y=fEr2(x)

	global gbl; #gbl=[gCa,gL,ECa,EL,parmT,parhT];

	m1r=fmh(0.,x,1,-1,gbl(5:9),0,0,0);   # steady-state activ. at Er	
	h1r=fmh(0.,x,1,-1,gbl(10:14),0,0,0);   # steady-state inactiv. at Er
	m3hr=m1r^3*h1r;
	y=gbl(1)*m3hr*(x-gbl(3))+gbl(2)*(x-gbl(4));

endfunction
