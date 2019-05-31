# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#	Steady-state current balance equation for computing Er, the
#	resting potential of a neurone:
#	   IH+IL=0
#	The function fEr2H below is used by fsolve a solver routine of
#	systems of nonlinear equations.
#	The parameters are passed as global var.s to the function.
#	The fct. fmh() computes the time course of m(t), h(t) or m(t)^k*h(t)
#	   of a membrane current when t, V(t) and init value(s) are given;



function y=fEr2H(x)

	global gbl; #gbl=[gH,gL,EH,EL,parmH];

	mHr=fmh(0.,x,1,-1,gbl(5:11),0,0,0);   # steady-state activ. at Er	
	mH3r=mHr^3;
	y=gbl(1)*mH3r*(x-gbl(3))+gbl(2)*(x-gbl(4));

endfunction
