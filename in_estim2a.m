# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#   This procedure computes the estimates: Cm, gL and EL for neurones when 
#   hyperpolarizing current steps are applied to them. In this version (2a),
#   it is assumed that the leakage and the capacitive current completely
#   determine the time evolution of the membrane pot. during the 
#   hyperpolarizing current steps.
#   The estimation procedure is the lin. least square (LLSQ) method. 
#   It is assumed that a constant holding current Ihold is applied to keep 
#   the cell in steady-state at the membrane potential Eh. Then the 
#   equilibrium condition is 
#		gL*(Eh-EL)=Ihold.			(1)
#   During the hyperpolarizing step
#		Cm*dV/dt+gL*(V-EL)=Ihold+Istep		(2)
#   applies. Subtracting (1) from (2), we get the equations for the LLSQ
#   method, namely one eqn. at each sampling time point. dV/dt is computed 
#   (approximated) via the Chebyshev approximation of V(t).
#
#   The neuronal input data is called burst, and its structure is given below.
#   In this version (2a), a correction value for a possible d.c. shift of 
#   the input currents must be supplied, in order to obtain the unbiased value 
#   of Ihold. The input current can also be re-scaled by the scaling factor 
#   cifc (cf. Input below).
#
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
#   k0: the data column of a voltage trace selected (1<k0<K0+2);
#   mch: no. of coeff.s  in the Chebyshev approx.;
#   dc_shift: d.c. shift of the current signal (nA); This is usually the mean
#	value of a current trace with no stimulus and no holding current;
#   cifc: correction factor for scaling the input current ci2.
#
#		Output:
#   Cm: membrane capacitance (pF);
#   gL: leakage conductance (nS);
#   cih: holding current applied (pA);
#   Eh: holding potential (mV) at cih;
#   EL: rev. pot. (mV) of the leakage current IL;
#   Er: resting potential (mV), Er=EL;
#   ci2: the constant stimulus current (pA);
#   t: vector of the selected sampled time (ms) on [tb,te];
#   v: vector of the selected sampled voltage data V corresponding to t (mV);
#   vch: Chebyshev approx. vector of the voltage trace V (mV) on [tb,te];
#   dv: Chebyshev approx. vector of dV/dt (mV/ms) on [tb,te];
#   IL: leakage curr. (pA) on [tb,te].
#   In addition sigma and the max. absolute error errmax=max(abs(Y-X*[Cm,gL]'))
#   of the LLSQ estimation are displayed (cf. Octave's User's Guide).
#
#		Internal variables:
#   Y: right hand side of the lin. equation;
#   X: coeff. matrix for Cm and gL.
#
#   


function [Cm,gL,cih,Eh,EL,Er,ci2,t,v,vch,dv,IL]=in_estim2a(burst,\
	tsb,tse,tb,te,k0,mch,dc_shift,cifc)

	N0=rows(burst);		# no. of data point in a trace

# First checking consistency of the data:
	if (N0<50) error("In in_estim2a: too few data points!\n") endif
	if ((tsb>=tse)||(tsb<burst(1,1))||(tse>burst(N0,1)))
	   printf("tsb=%f\t tse=%f\n",tsb,tse);
	   error("In in_estim2a: tsb or tse incorrect!\n")
	endif
	if ((tb<tsb)||(te>tse))
	   printf("tsb=%f\t tse=%f\t tb=%f\t te=%f\n",tsb,tse,tb,te);
	   error("In in_estim2a: tb or te is incorrect!\n")
	endif
	k1=floor((columns(burst)-1)/2);	   # no. of col.s of voltage data
	k2=columns(burst)-2*k1-1;
	if (k2>0) error("In in_estim2a: data are incomplete!\n") endif
	if ((k0==1)||(k0>k1+1))
	   k0
	   error("In in_estim2a: k0 is incorrect!\n")
	endif
	if (mch<3)
	   mch
	   error("In in_estim2a: order of Chebyshev approx. is too low!\n")
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

# LLSQ procedure:
	Y=ci2*ones(N,1);
	X=zeros(N,2);
# The following expression contains the equlibrium condition at Er:
	X(:,1)=dv;
	X(:,2)=v-Eh;
	[P,sig,R]=ols(Y,X);
	sig
	errmax=max(abs(R))		# max. abs. error of the LLSQ estim.
# Computing the output var.s:
	Cm=P(1);
	gL=P(2);

# Other output:
	EL=Eh-cih/gL;
	Er=EL;
	IL=gL*(v-EL);

endfunction
