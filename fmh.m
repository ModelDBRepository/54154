# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#   This routine computes the m^p*h^q(V(t)) function values for a 
#   given (t,V) vector pair, where V=V(t). m is the activation, 
#   h the inactivation var. of a voltage-gated ion current.
#   This routine implements the local solution formula.
#
#		Input:
#   t: time vector;
#   V: vector of the membrane pot.s as fct. of t;
#   p: order of the activ. kinetics (integer 1<=p<=4);
#   m0: init. val. of the activ. var.;
#   parm: par. vector of the activ. (cf. fm_gen(..));
#   q: inactiv. indicator: q=1 inactiv. current; q=0 non-inactiv. curr.;
#   h0: init. val. of the inactiv. var.;
#   parh: par. vector of the inactiv.: (cf. fm_gen(..));
#
#		Output:
#   xmph: vector of the values m(V(t))^p*h(V(t)) or m(V(t))^p.
#
#		External functions:
#  fm_gen(): computes the specific v-gated currents: IT,IA,IH,INa and IK
#	     the properties of which are specified in their par.vectors;


function xmph=fmh(t,V,p,m0,parm,q,h0,parh)

	N=length(t);
	if (N!=length(V))
	   error("In fmh: V and t are incompatible!\n");
	endif

# Step-by-step extrapolation of the values of m and h:
	t0=t(1);
	for ii=1:N
	   dt=t(ii)-t0;
	   v1=V(ii);
	   m1=fm_gen(v1,dt,m0,parm);
	   if (q==1) 
	      h1=fm_gen(v1,dt,h0,parh);
	      xmph(ii)=m1^p*h1;
	      h0=h1;
	   else
	      xmph(ii)=m1^p;
	   endif
	   t0=t(ii);
	   m0=m1;
	endfor

endfunction


