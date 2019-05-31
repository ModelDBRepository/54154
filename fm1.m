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
#	ionic current. The activation properties are given in the
#	alpha-beta form.
#
#		Input:
#   t: time vector;
#   V: vector of the membrane pot.s as fct. of t;
#   m0: init. val. of activ.;
#   paralp: vector of the alpha(V) parameters;
#   parbet: vector of the beta(V) parameters;
#   
#		Output:
#   xm: vector of the values m(V(t)).
#   xmse: vector of the steady-state val.s; 
#   gamme: vector of the gamma values on V;
#

function [xm,xmse,gamme]=fm1(t,V,m0,paralp,parbet)

	N=length(t);
	if ((N!=length(V)))
	   error("In fm1: V and t are incompatible!\n");
	endif


# Recursive extrapolation of the values of m:
	t0=t(1);
	for ii=1:N
	   dt=t(ii)-t0;
	   v1=V(ii);
	   v1a=paralp(2)*v1+paralp(3);
	   if (abs(v1a)<1e-6)
	      alp=paralp(4)+paralp(1)/paralp(2);
	   else
	      alp=paralp(1)*v1a./(exp(v1a)-1)+paralp(4);
	   endif
	   v1a=parbet(2)*v1+parbet(3);
	   if (abs(v1a)<1e-6)
	      bet=parbet(4)+parbet(1)/parbet(2);
	   else
	      bet=parbet(1)*v1a./(exp(v1a)-1)+parbet(4);
	   endif
	   gam1=alp+bet;
	   mss=alp/gam1;
	   m1=mss-(mss-m0)*exp(-dt*gam1);
	   xm(ii)=m1;
	   xmse(ii)=mss;
	   gamme(ii)=gam1;
	   t0=t(ii);
	   m0=m1;
	endfor	   
    
	  

endfunction


