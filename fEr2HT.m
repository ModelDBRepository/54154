# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#	Steady-state current balance equation for computing Er, the
#	resting potential of a neurone:
#	   IH+IT+IL=0
#	The function fEr2HT below is used by fsolve a solver routine of
#	systems of nonlinear equations.
#	The parameters are passed as global var.s to the function.
#	The fct. fmh() computes the time course of m(t), h(t) or m(t)^k*h(t)
#	   of a membrane current when t, V(t) and init value(s) are given;



function y=fEr2HT(x)

	global gbl; #gbl=[gH,gCa,gL,EH,ECa,EL,parmH,parmT,parhT];

	mHr=fmh(0.,x,1,-1,gbl(7:13),0,0,0);   # steady-state activ. at Er	
	mH3r=mHr^3;
	m1r=fmh(0.,x,1,-1,gbl(14:18),0,0,0);	# steady-state activ. at Er	
	h1r=fmh(0.,x,1,-1,gbl(19:23),0,0,0);	# steady-state inactiv. at Er
	m3hr=m1r^3*h1r;
	y=gbl(1)*mH3r*(x-gbl(4))+gbl(2)*m3hr*(x-gbl(5))+gbl(3)*(x-gbl(6));

endfunction
