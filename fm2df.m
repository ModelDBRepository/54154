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
#	This procedure is applied when the activ. curve (m_inf) is a product
#	of a (steep) Boltzmann curve and a polynomial of in V, and 
#	gamma(V) or 1/gamma(V) is a polynomial in V. 
#	The last entry of the par. vect. gamx indicates which is the
#	case. If this entry is positive (+1) then the former, if it is 
#	negative (-1) then the latter.
#	gamma(V) is bounded from below: gamma(V)>=1.
#
#		Input:
#   t: time vector;
#   V: vector of the membrane pot.s as fct. of t;
#   m0: init. val. of activ.;
#   p_xm: coeff. vector of the polynom. describing m_inf.;
#   q1: `slope factor' of the Boltzmann curve of m_inf(V);
#   V0: Vhalf of the Boltzmann curve of m_inf(V);
#   gamx: coeff. vector of the polynomial for gamma(V) or 1/gamma(V),
#	  depending on its last entry;
#
#		Output:
#   xm: vector of the values m(V(t)).
#   xmse: vector of the steady-state val.s; 
#   gamme: vector of the values gamma(V(t));
#
#		External functions:
#   minf(): Boltzmann curve;


function [xm,xmse,gamme]=fm2df(t,V,m0,p_xm,q1,V0,gamx)

	N=length(t);
	if ((N!=length(V)))
	   error("In fm2df: V and t are incompatible!\n");
	endif

	ngam=length(gamx);
	iflag=gamx(ngam);
	gamx=gamx(1:ngam-1);

# Recursive extrapolation of the values of m:
	t0=t(1);
	for ii=1:N
	   dt=t(ii)-t0;
	   v1=V(ii);
	   mss=max(1.e-8,polyval(p_xm,v1)*minf(v1,q1,V0));
	   gam1=polyval(gamx,v1);
	   if (iflag>0)
	      gam1=max(1,gam1);
	   else
	      gam1=max(1,1/gam1);
	   endif
	   m1=mss-(mss-m0)*exp(-dt*gam1);
	   xm(ii)=m1;
	   xmse(ii)=mss;
	   gamme(ii)=gam1;		# N.B. output is gamma(V)!
	   t0=t(ii);
	   m0=m1;
	endfor	   

endfunction


