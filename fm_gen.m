# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#   This routine computes the actual value xm(t) of an activ., or inactiv. 
#   var. at time t=t0+dt from its value m0 at time t0, using the param. 
#   values of the steady-state curve and the kinetics. If m0<0, the routine
#   returns the steady-state value at V.
#
#		Input:
#  V: the membrane pot. at t0 (it can be the vector of the spatial 
#	distribution);
#  d_t: scalar time increment;
#  m0: value of the var. at t0 (it can be the vector of the spatial 
#	distribution); 
#  par: param. vector of the steady-state and kinetc properties of the current:
#	1. in the (m_inf,tau)-type description of the current: par(1)=Vhalf; 
#	par(2)=k; par(3)=tau_min; par(4:5)=coeffs. of tau(V) as expon. fct.;
#	2. in the (alpha,beta)-type description of the current:
#	par(1:3)= par.s of alpha(V); par(4:6)=par.s of beta(V);
#	par(7)=current type, i.e. 1: mK; 2: mNa; 3: hNa;
#	3. IH: par(1:2): Botzmann steady-state activ. curve; par(4:6):
#	two branches of time const.s switching at par(7).
#	4. Inactiv. of INa according to Mainen et al. (Neuron, 15:1427,1995).
#	par(1:6): alpha-beta description of the rate const. as for the activ.
#	of INa; par(7:8): Boltmann steady-state activ. curve.
#
#		Output: 
#  xm: the updated values of m(t) (it can be the vector of the spatial 
#	distribution); 




function xm=fm_gen(V,d_t,m0,par)

	pl=length(par);

	if (pl==5)			# T current or A current
# steady-state curve:
           va=(V+par(1))/par(2);
           m1=1./(1+exp(va));		# Boltzmann steady-state curve
	   if (par(4)==0)
	      r_tb=1/par(3);		# kinetics is v-independent
	   elseif (par(5)==0)
	      r_tb=1/(par(3)+par(4));	# kinetics is v-independent
	   else
              va=-par(5)*V;
              r_tb=1./(par(3)+par(4)*exp(va));		# 1/(time const.)
	   endif
	elseif  (pl==7)				# currents IK, INa, or IH
	   if (par(7)==1)			# m_K
	      va=V+par(2);
	      if (abs(va)>1.e-5)
		 am=par(1)*va./(1-exp(-par(3)*va));
	      else
		 am=par(1)/par(3);
	      endif
	      bm=par(4)*exp(-par(6)*(V+par(5)));
	      r_tb=am+bm;
	      m1=am./r_tb;
	   elseif (par(7)==2)			# m_Na
	      va=V+par(2);
	      if (abs(va)>1.e-5)
		 am=par(1)*va./(1-exp(-par(3)*va));
	      else
		 am=par(1)/par(3);
	      endif
	      va=V+par(5);
	      if (abs(va)>1.e-5)
		 bm=par(4)*va./(exp(par(6)*va)-1);
	      else
		 bm=par(4)/par(6);
	      endif
	      r_tb=am+bm;
	      m1=am./r_tb;
	   elseif (par(7)==3)			# h_Na
	      am=par(1)*exp(-par(3)*(V+par(2)));
	      bm=par(4)./(1+exp(-par(6)*(V+par(5))));
	      r_tb=am+bm;
	      m1=am./r_tb;
	   else					# m_H of IH
              va=(V+par(1))/par(2);
              m1=1./(1+exp(va));	# Boltzmann steady-state curve
	      if (V>par(7))
	         r_tb=par(5)*exp(par(6)*V);
	      else
	         r_tb=par(3)*exp(par(4)*V);
	      endif
	      r_tb=1./r_tb;		# 1/(time const.)
	   endif
	elseif (pl==8)			# Mainen et al.s' h_Na of INa
	   va=V+par(2);
	   if (abs(va)>1.e-5)
	      am=par(1)*va./(1-exp(-par(3)*va));
	   else
	      am=par(1)/par(3);
	   endif
	   va=V+par(5);
	   if (abs(va)>1.e-5)
	      bm=par(4)*va./(exp(par(6)*va)-1);
	   else
	      bm=par(4)/par(6);
	   endif
	   r_tb=am+bm;			# alpha-beta rate const.
           va=(V+par(7))/par(8);
           m1=1./(1+exp(va));		# Boltzmann steady-state curve
	else
	   error("In fm_gen: This par. vector is invalid!\n");
	endif


	if (m0<0)
	  xm=m1;	# steady-state val. is returned
	  return
	endif

	if (d_t==0)
	   xm=m0;
	else
           xm=m1-(m1-m0).*exp(-d_t*r_tb);
	endif


endfunction
