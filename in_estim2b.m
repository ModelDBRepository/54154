# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#   This procedure computes the estimates: Cm, gCa, gL and EL for neurones when 
#   hyperpolarizing current steps are applied to them. In this version (2b),
#   it is assumed that the low-threshold Ca++ current IT, the leakage 
#   current IL, and the capacitive current completely determine the time 
#   evolution of the membrane pot. during the  hyperpolarizing current steps.
#   The estimation procedure is the lin. least square (LLSQ) method. 
#   It is assumed that a constant holding current Ihold is applied to keep 
#   the cell in steady-state at the membrane potential Eh.
#   Then the equilibrium condition is 
#	gCa*mT0^3(Eh)*hT0(Eh)*(Eh-ECa)+gL*(Eh-EL)=Ihold	(1)
#   where mT0(Eh) is the steady-state value of the activation var. mT of IT 
#   and hT0(Eh) is that of the inactivation var. hT of IT at Eh. ECa is 
#   the rev. pot. of IT.  
#   During the hyperpolarizing step
#	Cm*dV/dt+gCa*mT(t)^3*hT(t)*(Eh-ECa)+gL*(V-EL)=Ihold+Istep	(2)
#   applies with mT(t) and hT(t) as the time evolution of mT and hT,
#   respectively. These time fct.s are computed using the local solution 
#   formula (extrapolation) in the fct. fmh(). 
#   Subtracting (1) from (2), we get the equations for the LLSQ method, 
#   namely one eqn. at each sampling time point. dV/dt is computed 
#   (approximated) via the Chebyshev approximation of V(t). The resting 
#   potential Er is computed as the solution to (1) with Ihold=0, once gCa, 
#   gL and EL are known.
#
#   The neuronal input data is called burst, and its structure is given below.
#   In this version (2b), a correction value for a possible d.c. shift of 
#   the input currents must be supplied, in order to obtain the unbiased value 
#   of Ihold. The input current can also be re-scaled by the scaling factor 
#   cifc (cf. Input below).
#
#		Input:
#   burst: data matrix whose 1st column is the sampled time; the subsequent
#	   K0 columns are the sampled voltage at different hyperpol. stimuli;
#	   and the last K0 columns are the sampled stimuli corresponding to
#	   the voltage traces, i.e. columns i and i+K0 correspond to each
#	   other (i>1). The no. of rows is the length of a data trace. The
#	   voltage values are given in mV, those of the current in nA;
#   tsb: start time of the stimulus; all time data are given in ms;
#   tse: end of the stimulus;
#   tb: start time of the data segment used in the estimation, tb>tsb;
#   te: end time of the data segment used in the estimation, te<tse;
#   k0: data column of voltages selected (1<k<K0+2);
#   mch: no. of coeff.s in the Chebyshev approx.;
#   dc_shift: d.c. shift of the current signal (nA); This is usually the mean
#	value of a current trace with no stimulus and no holding current;
#   cifc: correction factor for scaling the input current ci2;
#   parmT: par. vector of the activ. of IT;
#   parhT: par. vector of the inactiv. of IT;
#   ECa: reversal pot. for Ca (mV);
#
#		Output:
#   Cm: membrane capacitance (pF);
#   gCa: maximal T-conductance (nS);
#   gL: leakage conductance (nS);
#   cih: holding current applied (pA);
#   Eh: holding potential (mV) at cih;
#   EL: rev. pot. (mV) of IL , the leakage current;
#   Er: resting potential (mV), Er=EL;
#   ci2: the average constant stimulus current (pA);
#   t: vector of the selected sampled time (ms);
#   v: vector of the selected sampled voltage data corresponding to t (mV);
#   vch: Chebyshev approx. vector of the voltage trace (mV);
#   dv: Chebyshev approx. vector of dV/dt (mV/ms);
#   IL: leakage curr. (pA) on [tb,te];
#   IT: T-curr. (pA) on [tb,te];
#   mT: m(t) of IT on [tb,te];
#   hT: h(t) of IT on [tb,te];
#   In addition sigma and the max. absolute error of the LLSQ estimation 
#	are displayed (cf. Octave's User's Guide).
#
#		Internal variables:
#   m1r: value of the activ. var. of IT at Eh;
#   h1r: value of the inactiv. var. of IT at Eh;
#   m3hr: m1r^3*h1r;
#   m3h: vector of the m^3*h values at the sampled (v,t) pairs;
#   Y: right hand side of the lin. equation;
#   X: coeff. matrix for Cm, gCa and gL;
#
#		External functions and procedures:
#   fmh(): computes the time course of m(t), h(t) or m(t)^k*h(t)
#	   of a membrane current when t, V(t) and init value(s) are given;
#   fEr2(): current balance equation used in fsolve() to compute Er (cf.
#	   Octave's User's Guide).
#



function [Cm,gCa,gL,cih,Eh,EL,Er,ci2,t,v,vch,dv,IL,IT,mT,hT]=\
	in_estim2b(burst,tsb,tse,tb,te,k0,mch,dc_shift,cifc,parmT,parhT,ECa)

	N0=rows(burst);		# no. of data point in a trace

# First checking consistency of the data:
	if (N0<50) error("In in_estim2b: too few data points!\n") endif
	if ((tsb>=tse)||(tsb<burst(1,1))||(tse>burst(N0,1)))
	   printf("tsb=%f\t tse=%f\n",tsb,tse);
	   error("In in_estim2b: tsb or tse incorrect!\n")
	endif
	if ((tb<tsb)||(te>tse))
	   printf("tsb=%f\t tse=%f\t tb=%f\t te=%f\n",tsb,tse,tb,te);
	   error("In in_estim2b: tb or te is incorrect!\n")
	endif
	k1=floor((columns(burst)-1)/2);	   # no. of col.s of voltage data
	k2=columns(burst)-2*k1-1;
	if (k2>0) error("In in_estim2b: data are incomplete!\n") endif
	if ((k0==1)||(k0>k1+1))
	   k0
	   error("In in_estim2b: k0 is incorrect!\n")
	endif
	if (mch<3)
	   mch
	   error("In in_estim2b: order of Chebyshev approx. is too low!\n")
	endif

# Compute Eh (holding pot.) as average pre-stim. potential:
	for n0=1:N0
	  if (tsb<=burst(n0,1)) break; endif  # start time of stimulus
	endfor
	Eh=mean(burst(10:n0-20,k0));	      # avoid transients

# Select data segment for the LLSQ, and find the end of the stimulus:
	for n1=n0:N0
	  if (tb<=burst(n1,1)) break; endif  # start time of selected data
	endfor
	for n2=n1:N0
	  if (te<=burst(n2,1)) break; endif  # end time point of selected data
	endfor
	for n3=n2:N0
	  if (tse<=burst(n3,1)) break; endif  # end of stimulus
	endfor
	N=n2-n1+1;		# no. of data points in the selected segment
	t=burst(n1:n2,1);
	v=burst(n1:n2,k0);

# Compute the stim. curr. ci2 as average of the sampled values 
# in [tsb,tse], and the pre-stimulus holding curr. cih
# avoiding transients:
	ci2=mean(burst(n0+10:n3-10,k0+k1));
	cih=mean(burst(10:n0-20,k0+k1));
	ci2=1000*(ci2-cih);			# in pA
# Correct the pre-stimulus holding current cih for dc shift:
	cih=1000*(cih-dc_shift);		# in pA
# Re-scale the input current if necessary:
	if ((cifc !=1)&&(cifc !=0))
	   ci2=cifc*ci2;
	   cih=cifc*cih;
	endif

# Compute the Chebyshev approx. of the voltage data and that of their
# time derivative dV/dt:
	ck=cheb_linip(t,v,mch);		# Chebyshev coeff.s 
	vch=chebev_vect(t(1),t(N),ck,t)';  
	dv=df_ch_vect(t,v,mch)';

# Compute the activ. and inactiv. var.s of IT:
	m1r=fmh(0.,Eh,1,-1,parmT,0,0,0);   # steady-state activ. at Eh	
	h1r=fmh(0.,Eh,1,-1,parhT,0,0,0);   # steady-state inactiv. at Eh
	m3hr=m1r^3*h1r;
	m3h=fmh(t,v,3,m1r,parmT,1,h1r,parhT);   # for all (t,v);
	hT=fmh(t,v,1,h1r,parhT,0,0,0);		# h(t) of IT for all (t,v)
	mT=(m3h./hT).^(1/3);			# m(t) of IT for all (t,v)

# LLSQ procedure:
	Y=ci2*ones(N,1);
	X=zeros(N,3);
# The following expression contains the equlibrium condition at Eh:
	X(:,1)=dv;
	X(:,2)=m3h.*(v-ECa)-m3hr*(Eh-ECa);
	X(:,3)=v-Eh;
	[P,sig,R]=ols(Y,X);
	sig
	errmax=max(abs(R))		# max. abs. error of the LLSQ estim.
# Computing the output var.s:
	Cm=P(1);
	gCa=P(2);
	gL=P(3);

# Other output:
	EL=Eh+(gCa*m3hr*(Eh-ECa)-cih)/gL;
	IL=gL*(v-EL);
	IT=gCa*m3h.*(v-ECa);

# Resting potential Er:
	gbl=[gCa,gL,ECa,EL,parmT,parhT];
	global gbl;
	[Er,info]=fsolve("fEr2",Eh);
	if (info!=1) perror("fsolve",info); endif


endfunction
