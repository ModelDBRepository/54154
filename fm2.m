# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This routine computes the m(t) (or h(t)) function by extrapolation
#	for a given (t,V) vector pair, where V=V(t).
#	m is the activation, or inactivation var. of a voltage-gated
#	ionic current.
#	The steady-state activation is a Boltzmann curve, and 
#	gamma(V)=1/tau(V) is a polynomial in V which is bounded from below
#	as indicated.
#
#		Input:
#   t: time vector;
#   V: vector of the membrane pot.s as fct. of t;
#   m0: init. val. of activ.;
#   q: "slope" factor of the s-s activ. (Boltzmann) curve of m;
#   V0: "half value" of the s-s activ. (Boltzmann) curve of m;
#   gamx: vector of coeffs of a polynomial that describes gamma(V);
#
#		Output:
#   xm: vector of the values m(V(t)).
#   xmse: vector of the steady-state val.s (points of the Boltzmann curve); 
#   gamme: vector of the gamma values on V;
#



function [xm,xmse,gamme]=fm2(t,V,m0,q,V0,gamx)

	N=length(t);
	if ((N!=length(V)))
	   error("In fm2: V and t are incompatible!\n");
	endif

# Recursive extrapolation of the values of m:
	t0=t(1);
	for ii=1:N
	   dt=t(ii)-t0;
	   v1=V(ii);
	   mss=1/(1+exp(q*(v1+V0)));
	   gam1=max(polyval(gamx,v1),0.1);
	   m1=mss-(mss-m0)*exp(-dt*gam1);
	   xm(ii)=m1;
	   xmse(ii)=mss;
	   gamme(ii)=gam1;
	   t0=t(ii);
	   m0=m1;
	endfor	   
    

endfunction


