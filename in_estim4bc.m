# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	The main purpose of this procedure is to compute ci0=IK+INa, where
#	INa and IK are the Na+ and K+ currents, respectively, shaping the
#	action potential (AP). In addition, a number of values and time 
#	functions as listed in the `Output:' below are also calculated.
#	It is assumed that the required properties of A-type K+ current IA, 
#	the low-threshold Ca++ current IT, the hyperpolarisation-activated
#	current IH, and the leakage current IL, as well as Cm, are known. 
#	The stimulus evokes at least one AP. In the computations, only the 
#	first AP is used and the time interval [tb,te] is set such that 
#	fully to include the first AP but to exclude subsequent APs.
#	The procedure computes the time course of the currents IA, IT, IH 
#	by applying the local solution formula (extrapolation from one 
#	sampled time point to the next) to solve the corresponding kinetic 
#	eqn.s. In this procedure, the activation and kinetic param.s as fct.s 
#	of V  obtained from v-clamp data are used. The current IA or IT can
#	be omitted by setting gA0=0 or gCa0=0. However, IH and IL are assumed
#	to be present unconditionally.
#	The derivative dV/dt is computed via Chebyshev approximation of V(t). 
#	Because of the fast changes of V(t) during the initial phase of the AP, 
#	the time interval is subdivided into a pre-AP segment and 
#	the AP-segment. dV/dt is computed separately for both segments. 
#	To ensure continuity of dV/dt between the segments, a third segment
#	overlapping each of the other segments is used in which dV/dt is 
#	again computed. The values of dV/dt in the overlapping segment serve 
#	to obtain an improved dV/dt by simple smoothing (averaging).
#	Note on choice of mch: Increase mch until the changes of the
#	current ci0=INa+IK between consecutive mch values become sufficiently 
#	small. Usually mch>20. (See also User's Guide.)
#	The injected current is corrected for dc shift and can be re-scaled.	
#
#		Input:
#   burst: data matrix whose 1st column is the sampled time; the subsequent
#	   K0 columns are the sampled time courses of the voltage at 
#	   different stimuli; and the last K0 columns are the sampled 
#	   stimuli corresponding to the voltage traces, i.e. columns i 
#	   and i+K0 correspond to each other (i>1). The no. of rows is 
#	   the length of a data trace. The voltage values are given in mV, 
#	   those of the current in nA;
#   tsb: start time of the stimulus; all time data are given in ms;
#   tse: end of the stimulus;
#   tb: start time of the data segment used in the estimation, tb>tsb;
#   te: end time of the data segment used in the estimation, te<tse;
#   k0: data column of voltage traces selected (1<k<K0+2);
#   mch: no. of coeffs. in the Chebyshev approx. used for V(t);
#   parmA: par. vector of the activ. of IA;
#   parhA: par. vector of the inactiv. of IA;
#   parmT: par. vector of the activ. of IT;
#   parhT: par. vector of the inactiv. of IT;
#   parmH: par. vector of the activ. of IH;
#   parhNa: par. vector of the inactiv. of INa;
#   EK: reversal pot. for K (mV);
#   ECa: reversal pot. for Ca (mV);
#   EH: rev. pot for the IH (mV);
#   EL: rev. pot for the leakage (mV);
#   Cm0: membrane capacitance (pF);
#   gA0: max. conductance of IA (nS);
#   gCa0: max. conductance of IT (nS);
#   gH0: max. conductance of IH (nS);
#   gL0: leakage conductance  (nS);
#   cifc: current scaling factor;
#
#		Output:
#   Er: the (averaged) resting potential (mV);
#   ci2: the averaged constant stimulus current (pA);
#   t: vector sampling times on [tb,te];
#   v: vector of sampled voltage data (mV) on [tb,te];
#   vch: Chebyshev approx. vector of the voltage trace (mV);
#   dv: Chebyshev approx. vector of dV/dt (mV/ms) on [tb,te];
#   ci0: INa(t)+IK(t) on [tb,te];
#   IA: IA(t) on [tb,te];
#   IT: IT(t) on [tb,te];
#   hNa: inactivation hNa(t) of INa on [tb,te];
#   ts1: time vector on [tsb,tb];
#   vs1: vector of smoothed sampled membrane pot. on [tsb,tb];
#   dvs1: time derivative of vs1 on [tsb,tb];
#   ci0s1: INa(t)+IK(t) on [tsb,tb];
#   IAs1: IA(t) on [tsb,tb];
#   ITs1: IA(t) on [tsb,tb];
#   hNas1: hNa(t) on [tsb,tb];
#
#		Internal variables:
#   ts0: `overlapping' time segment;
#   vs0: voltage values over ts0;
#   m1r: steady-state value of the activ. var. of IA or IT at Er;
#   h1r: steady-state value of the inactiv. var. of IA or IT at Er;
#   m4hr: m1r^4*h1r of IA;
#   m4hs: actual value of mA^4*hA(t) of IA on the interval [tsb,te];
#   IAs: actual value of IA on the interval on [tsb,te];
#   m3hr: m1r^3*h1r of IT;
#   m3hs: actual value of m^3*h(t) of IT on the interval [tsb,te];
#   ITs: actual value of IT on the interval on [tsb,te];
#   mHr: steady-state value of the activ. var. of IH at Er;
#   mH3s: actual value of mH^3(t) of IH on [tsb,te];
#   hNas: actual value of hNa(t) of INa on [tsb,te];
#   ci0s: INa(t)+IK(t) on [tsb,te];
#
#		External functions and procedures:
#   fmh(): computing the local solutions for activation-inactivation;




function [Er,ci2,t,v,vch,dv,ci0,IA,IT,hNa,ts1,vs1,dvs1,ci0s1,IAs1,ITs1,\
	hNas1]=in_estim4bc(burst,tsb,tse,tb,te,k0,mch,parmA,parhA,\
	parmT,parhT,parmH,parhNa,EK,ECa,EH,EL,Cm0,gA0,gCa0,gH0,gL0,cifc)


	N0=rows(burst);		# no. of data points in a trace

# Checking consistency of the data:
	if (N0<50) error("In in_estim4bc: too few data points!\n") endif
	if ((tsb>=tse)||(tsb<burst(1,1))||(tse>burst(N0,1)))
	   printf("tsb=%f tdat(1)=%f tse=%f tdat(last)=%f\n",\
		tsb,burst(1,1),tse,burst(N0,1));
	   error("In in_estim4bc: tsb or tse incorrect!\n")
	endif
	if ((tb<tsb)||(te>tse))
	   printf("tsb=%f\t tse=%f\t tb=%f\t te=%f\n",tsb,tse,tb,te);
	   error("In in_estim4bc: tb or te is incorrect!\n")
	endif
	k1=floor((columns(burst)-1)/2);	   # no. of col.s of voltage data
	k2=columns(burst)-2*k1-1;
	if (k2>0) error("In in_estim4bc: data are incomplete!\n") endif
	if ((k0==1)||(k0>k1+1))
	   k0
	   error("In in_estim4bc: k0 is incorrect!\n")
	endif
	if (mch<5)
	   mch
	   error("In in_estim4bc: order of Chebyshev approx. is < 4!\n")
	endif

# Compute Er as average pre-stim. potential:
	for n0=1:N0
	  if (tsb<=burst(n0,1)) break; endif  # start time of stimulus
	endfor
	Er=sum(burst(10:n0-10,k0))/(n0-19);   # avoid transients

# Select data segment, and find the end of the stimulus:
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
	ck=cheb_linip(t,v,mch);			# Cheb. coeffs of v
	vch=chebev_vect(t(1),t(N),ck,t)';	# Cheb. approx. of v
	dv0=df_ch_vect(t,v,mch)';		# temporary dV/dt

	Ns1=n1-n0+1;			# no. of data points in [tsb,tb]
	ts1=burst(n0:n1,1);
	vs1=burst(n0:n1,k0);

	mch2=5;				# low order Cheb. approx. suffices
	dvs1=df_ch_vect(ts1,vs1,mch2)';

# Find the `break point' where V starts increasing steeply on [tb,te]
	d0v1=dvs1(Ns1)*(t(3)-t(1));
	for nc0=2:N-1
	   d0v2=v(nc0+1)-v(nc0-1);
	   if (d0v2>20*d0v1) break; endif
	endfor

	if (nc0>=n2) error("In in_estim4bc: No action pot.!\n"); endif
	ns0=11;				  # length of segment used in vs1
	vs0=[vs1(Ns1-ns0+1:Ns1);v(2:nc0+1)];
	ts0=[ts1(Ns1-ns0+1:Ns1);t(2:nc0+1)];
	mch2=5;				# low order Cheb. approx. suffices
	dvs0=df_ch_vect(ts0,vs0,mch2)';

# Correct (smooth) dv0 by using dvs0
	dv=dv0;
	for ii=nc0:-1:1
	   if (dv0(ii)>0)
	      dv(ii)=0.5*(dv0(ii)+dvs0(ns0+ii-1));    # use average value
	   else
	      break;
	   endif
	endfor
	dv1=dvs1(Ns1);			# boundary (overlap) condition
	dv(1)=dv1;
	dv2=dv(ii+1);
	for jj=2:ii
	   dv1=((ii+1-jj)*dv1+dv2)/(ii+2-jj);	# new dv values by interpol.
	   dv(jj)=dv1;
	endfor

# Check and correct (inital part of) dv for monotonicity:
	inmx=min(find(1-sign(max(dv)-dv)));
	dv1=dv(1);
	for ii=2:inmx
	   ii1=ii-1;
	   dv2=dv(ii);
	   if ((dv2>dv1)&&(dv2>dv(ii1))) break; endif
	endfor
	nc0=ii;				# re-define nc0
	if (nc0>2)
	   t1=t(1);
	   for ii=nc0:-1:2
	      ii1=ii-1;
	      d0v2=(dv(ii)-dv1)/(t(ii)-t1);
	      dv(ii1)=d0v2*(t(ii1)-t1)+dv1;
	   endfor
	endif
	
# Compute the stim. curr. (in pA) as average of the sampled values 
# in [tsb,tse], and avoid transients:
	ci2=1000*mean(burst(n0+10:n3-10,k0+k1));
# Correct for dc shift in the current, and re-scale:
	ci20=1000*mean(burst(10:n0-20,k0+k1));	# dc shift
	if ((cifc !=1)&&(cifc !=0))
	   ci2=cifc*(ci2-ci20);
	endif

# Total time and voltage:
	ts=[ts1;t(2:N)];
	vs=[vs1;v(2:N)];
	dvs=[dvs1;dv(2:N)];

# Compute the activ. and inactiv. var.s of IA from v-clamp data if present:
	if (gA0>0)
	   m1r=fmh(0.,Er,1,-1,parmA,0,0,0);   # steady-state activ. at Er	
	   h1r=fmh(0.,Er,1,-1,parhA,0,0,0);   # steady-state inactiv. at Er
	   m4hr=m1r^4*h1r;
	   m4hs=fmh(ts,vs,4,m1r,parmA,1,h1r,parhA);   # for all (ts,vs);
# Compute IAs :
	   IAs=gA0*m4hs.*(vs-EK);
	   IAs1=IAs(1:Ns1);		       # restrict to (ts1,vs1)
	   IA=IAs(Ns1:Ns1+N-1);		       # restrict to (t,v)
	else 				       # no IA
	   IAs=0;
	   IAs1=0;
	   IA=0;
	endif

# Compute the activ. and inactiv. var.s of IT from v-clamp data if present:
	if (gCa0>0)
	   m1r=fmh(0.,Er,1,-1,parmT,0,0,0);   # steady-state activ. at Er	
	   h1r=fmh(0.,Er,1,-1,parhT,0,0,0);   # steady-state inactiv. at Er
	   m3hr=m1r^3*h1r;
	   m3hs=fmh(ts,vs,3,m1r,parmT,1,h1r,parhT);   # for all (ts,vs);
# Compute ITs :
	   ITs=gCa0*m3hs.*(vs-ECa);
	   ITs1=ITs(1:Ns1);		       # restrict to (ts1,vs1)
	   IT=ITs(Ns1:Ns1+N-1);		       # restrict to (t,v)
	else 				       # no IT
	   ITs=0;
	   ITs1=0;
	   IT=0;
	endif

# Compute the activ. var of IH from v-clamp data:
	mHr=fmh(0.,Er,1,-1,parmH,0,0,0);     # steady-state activ. at Er	
	mH3r=mHr^3;
	mH3s=fmh(ts,vs,3,mHr,parmH,0,0,0);   # for all (ts,vs);

# Compute the inactiv. var. of INa from v-clamp data:
	hNar=fmh(0.,Er,1,-1,parhNa,0,0,0);     # steady-state activ. at Er	
	hNas=fmh(ts,vs,1,hNar,parhNa,0,0,0);   # for all (ts,vs)
	hNas1=hNas(1:Ns1);		       # restrict to (ts1,vs1)
	hNa=hNas(Ns1:Ns1+N-1);		       # restrict to (t,v)

# Compute ci0=INa+IK
	ci0s=ci2-Cm0*dvs-gL0*(vs-EL)-gH0*mH3s.*(vs-EH)-IAs-ITs;
	ci0s1=ci0s(1:Ns1);		# restrict to (ts1,vs1)
	ci0=ci0s(Ns1:Ns1+N-1);		# restrict to (t,v)


endfunction   
