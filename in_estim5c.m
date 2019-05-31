# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This procedure computes the inital estimates of gNa, and the par.s of
#	alpha_Na(V) and beta_Na(V) in the form, shown below.
#	It uses the same order of Chebyshev approx. as during the computaion
#	of dV/dt in in_estim4bc(). 
#   Notation:  tb=t(1)>tsb and te=t(length(t))<tse 
#		(tsb: stim. start, tse: stim. end)
#
#
#		Input:
#   t: time vector on [tb,te];
#   v: sampled V(t) on [tb,te];
#   ci0: ci0=IK+INa on [tb,te];
#   mch: no. of polynom. terms in the Chebyshev approx. of V;
#   ts1: sampled time on [tsb,tb]; 
#   vs1: sampled V(t) on [tsb,tb];
#   Er: resting pot.;
#   hNa: hNa(t) on [tb,te];
#   hNas1: hNa(t) on [tsb,tb];
#   ENa: Na rev.pot.;
#   gammn: pre-defined minimal value of gamma_Na(V)=1/taum_Na(V) (activ.rate);
#
#		Output:
#   ka:  minimal index in vector t where ci0<0, index of tka;
#   icmn: index of tcmn where ci0 is minimal (has a negative peak);
#   gNa: estimated  max. conductance of INa;
#   alpNa0: vector of non-negative alpha_Na values on [tka,tcmn];
#   valp: corresponding V vector;
#   betNa0: vector of non-negative beta_Na values on [tka,tcmn];
#   vbet: corresponding V vector;
#   lambda: parameter in the solution to alpha_Na and beta_Na;
#   alpNae: alpha_Na(V) estimated as ap*V1/(exp(V1)-1) where V1=alp1*V+alp2 
#	    on V=v;
#   p_alpNa: param. vector [ap,alp1,alp2];
#   betNae: beta_Na(V) estimated as bp*V1/(exp(V1)-1) where V1=bet1*V+bet2 
#	    on V=v;
#   p_betNa: param. vector [bp,bet1,bet2];
#   mNaes1: mNa(t) on [tsb,tb];
#   mNass1: mNa_inf(V(t)) on [tsb,tb];
#   gamNas1: gamma_Na(V(t)) on [tsb,tb];
#   mNae: mNa(t) on [tb,te];
#   mNas: mNa_inf(V(t)) on [tb,te];
#   gamNa: gamma_Na(V(t)) on [tb,te];
#   INas1: reconstructed INa(t) on [tsb,tb];
#   INa: reconstructed INa(t) on [tb,te];
#   
#
#		Internal variables:
#   ci0mn: minimum of ci0(t);
#   tcmn: tcmn=t(icmn), end point of the estimation interval [tka,tcmn];
#   tka: tka=t(ka) starting point of the estimation interval [tka,tcmn];
#   N1: length of the selected sub-interval [tka,tcmn], N1=icmn-ka+1;
#   mch2: the same as mch;
#   a: a(t)=[INa/(hNa*(v-ENa))]^(1/3)=gNa^(1/3)*mNa(t) on [tka,tcmn];
#   c_a: vector of Cheb. coeffs. of a;
#   ach: Chebyshev approximate of a(t) on [tka,tcmn];
#   av_err_a_cheb: `average' error of the Cheb approx. of a(t);
#   Ha: coeff. matrix using c_a;
#   D: derivative matrix of the Cheb. approx.;
#   g0Na: g0Na=gNa^(1/3);
#   c_alpNa0: vector of Cheb. coeffs. of alpha_Na0;
#   c_betNa0: vector of Cheb. coeffs. of beta_Na0;
#   G: coeff. matrix of polynom. lin. lsq.;
#   p_alp: coeff vector [alp1,alp2];
#   p_bet: coeff vector [bet1,bet2];
#   sig, res: sigma val. and residual vector of the lin. least square proc.;
#   m0: steady-state value of mNa at Er;
#   m00: dummy var.;
#
#
#		External procedures and functions:
#  cheb_linip(): computes the Cheb. coeffs.; 
#  chebev_vect(): computes fct values from the Cheb. coeffs.;
#  hgv3(): creates coeff. mtx from the Cheb coeffs.; 
#  d_cch(): derivative matrix of the Cheb. coeffs.;
#  f_eq1(): nonlin.eqn. given as 1+b*x=exp(x) used in fsolve(); 
#  fm1(): computes the actual values of the (in)activ. over a time vector 
#	   t when alpha(V) and beta(V) are a*x/(exp(x)-1), where x=b1*V+b2,
#	   as well as the time course of m_inf(V(t)) and gamma(V(t))=1/tau;





function [ka,icmn,gNa,alpNa0,valp,betNa0,vbet,lambda,alpNae,p_alpNa,\
	betNae,p_betNa,mNaes1,mNass1,gamNas1,mNae,mNas,gamNa,INas1,INa]=\
	in_estim5c(t,v,ci0,mch,ts1,vs1,Er,hNa,hNas1,ENa,gammn)


	N0=length(t);		# no. of data points on [tb,te]
# Check consistency of the data:
	if ((length(v) !=N0)||(length(ci0)!=N0)||(length(hNa)!=N0))
	   printf("data lengths: lt=%3d lv=%3d lci0=%3d lhNa=%3d\n",\
		N0,length(v),length(ci0),length(hNa));
	   error("In in_estim5c: Data are inconsistent\n");
	endif

# Find the negative peak of ci0:
	ci0mn=min(ci0);
	icmn=find(1+sign(ci0mn-ci0));
	tcmn=t(icmn);			# end point of the estim. interval

	for k=icmn:-1:1
	   if (ci0(k)>=0) break; endif
	endfor
	ka=k+1;
	tka=t(ka);			# starting point of the estim. interval
	N1=icmn-ka+1;			# no. of points in the estim. interv.
	mch2=mch;			# the same as used for V;
	a=ci0(ka:icmn)./(hNa(ka:icmn).*(v(ka:icmn)-ENa));	# gNa*mNa^3
	a=a.^(1/3);
	c_a=cheb_linip(t(ka:icmn),a, mch2);	# Cheb.coeff.vect of a
	ach=chebev_vect(tka,tcmn,c_a,t(ka:icmn))';
	av_err_a_cheb=max(abs(a-ach))*length(a)/sum(a)	# average error

	Ha=hgv3(mch2-1,c_a);			
	Ha=Ha(1:mch2,:);
	D=d_cch(tka,tcmn,mch2);			# derivative matrix
	dca=D*c_a;

	g0Na=1.15*max(eig(Ha));			# 1st estimate of gNa^(1/3)
	c_alpNa0=(g0Na*eye(mch2)-Ha)\dca;
	c_betNa0=Ha\dca;
	alpNa0=chebev_vect(tka,tcmn,c_alpNa0,t(ka:icmn))';
	betNa0=chebev_vect(tka,tcmn,c_betNa0,t(ka:icmn))';


# Use only positive values for alpNa0 and betNa0 (select the corresponding
# values of v accordingly):
	ka2=0; kb2=0;
	for k=1:N1
	   al0=alpNa0(k);
	   if (al0>0)
	      ka2++;
	      alpNa0(ka2)=al0;
	      valp(ka2)=v(ka+ka2-1);
	   endif
	   be0=betNa0(k);
	   if (be0>0)
	      kb2++;
	      betNa0(kb2)=be0;
	      vbet(kb2)=v(ka+kb2-1);
	   endif
	endfor
	alpNa0=alpNa0(1:ka2);
	betNa0=betNa0(1:kb2);

	al0=0;
	be0=0;

# Estimate par. values of alpha_Na and beta_Na;
# Coeff. matrix of the lin.sq. est. for alpha_Na:
	G(:,2)=ones(size(valp));
	G(:,1)=valp;
# First solve a nonlin.eqn. at each value of alpNa0:
	global ap ri;
	sq_err0=1e10;
	dcf0=0.5*abs(min(alpNa0(2:ka2)-alpNa0(1:ka2-1)));
	cf0mx=2*(max(alpNa0)+dcf0-al0)+1e-10;
	cf0=0.7*(min(alpNa0)-al0);
   while(cf0<cf0mx)
	ap=cf0;
	for ii=1:ka2
	   ri=alpNa0(ii)-al0;
	   if (ap==ri)
	      y0=0;
	      info=1;
	   elseif (ap>ri)
	      [y0,info]=fsolve("f_eq1",[10]);
	   else
	      [y0,info]=fsolve("f_eq1",[-ri/ap]);
	   endif
	   if (info!=1)
	      perror("fsolve",info)
	      printf("Error during estim. of alpNa0 at cf0=%5.2f\n",cf0);
	      break;
	   endif
	   ya(ii)=y0;
	endfor
	[p_alp,sig,res]=ols(ya,G);
# Compute the square error of the estimation
	va1=polyval(p_alp,valp);
	alpNa1=ap*va1./(exp(va1)-1);
	sq_err1=sumsq(alpNa0-al0-alpNa1);
	if (sq_err1<sq_err0)
	   p_alp0=p_alp;
	   sig_alp=sig;
	   res0=max(abs(res));
	   ap0=ap;
	   sq_err0=sq_err1;
	   cf0a=cf0;
	endif
	cf0=cf0+dcf0;
   endwhile
	printf("al0=%6.3f ap0=%6.3f alp1=%5.2f alp2=%5.2f\n",al0,ap0,p_alp0);
	printf("sig_alp=%13.6g res=%13.6g\n",sig_alp,res0);

# Now solve a nonlin. eqn. at each value of betNa0:
# Re-define coeff. matrix of the lin.sq. est. for beta_Na:
	clear G ya;
	G(:,2)=ones(size(vbet));
	G(:,1)=vbet;
	sq_err0=1e10;
	dcf0=0.5*abs(min(betNa0(2:kb2)-betNa0(1:kb2-1)));
	cf0mx=2*(max(betNa0)+dcf0-be0)+1e-10;
	cf0=0.7*(min(betNa0)-be0);
   while(cf0<cf0mx)
	ap=cf0;
	for ii=1:kb2
	   ri=betNa0(ii)-be0;
	   if (ap==ri)
	      y0=0;
	      info=1;
	   elseif (ap>ri)
	      [y0,info]=fsolve("f_eq1",[10]);
	   else
	      [y0,info]=fsolve("f_eq1",[-ri/ap]);
	   endif
	   if (info!=1)
	      perror("fsolve",info)
	      printf("Error during estim. of betNa0 at cf0=%5.2f\n",cf0);
	      break;
	   endif
	   ya(ii)=y0;
	endfor
	[p_bet,sig,res]=ols(ya,G);
# Compute the square error of the estimation
	va1=polyval(p_bet,vbet);
	betNa1=ap*va1./(exp(va1)-1);
	sq_err1=sumsq(betNa0-be0-betNa1);
	if (sq_err1<sq_err0)
	   p_bet0=p_bet;
	   sig_bet=sig;
	   res0=max(abs(res));
	   bp0=ap;
	   sq_err0=sq_err1;
	   cf0b=cf0;
	endif
	cf0=cf0+dcf0;
   endwhile
	printf("be0=%6.3f bp0=%6.3f bet1=%5.2f bet2=%5.2f\n",be0,bp0,p_bet0);
	printf("sig_bet=%13.6g res=%13.6g\n",sig_bet,res0);
 

# Compute the estimated alpha_Na0(V) and beta_Na0(V) for V=v:
	va2=polyval(p_alp0,v);
	alpNa0e=ap0*va2./(exp(va2)-1)+al0;
	va2=polyval(p_bet0,v);
	betNa0e=bp0*va2./(exp(va2)-1)+be0;

# Set lambda=10:
	lambda=10;
	alpNae=(1+lambda)*alpNa0e;
	betNae=lambda*betNa0e;
# Compute new lambda such that gam0=gammn:
	gam0=min(alpNae+betNae);
	igam0=min(find(1+sign(gam0-(alpNae+betNae))));
# Re-define alphaNae and betaNae with the new lambda:
	lambda=lambda+(gammn-gam0)/(alpNa0e(igam0)+betNa0e(igam0));
	alpNae=(1+lambda)*alpNa0e;
	betNae=lambda*betNa0e;
# Define par. vectors of alpNa_(V), beta_Na(V) for computing mNa(t):
	p_alp2=[(1+lambda)*ap0;p_alp0;(1+lambda)*al0];
	p_bet2=[lambda*bp0;p_bet0;lambda*be0];
	p_alpNa=p_alp2(1:3);	# usual par. set of alpha_Na(V)
	p_betNa=p_bet2(1:3);	# usual par. set of beta_Na(V)

# Compute mNa(t), mNa_inf(V(t)), and gamma_Na(V(t)) during ts1 and t:	
	[m00,m0,gamNa0]=fm1(0,Er,0,p_alp2,p_bet2);
	[mNaes1,mNass1,gamNas1]=fm1(ts1,vs1,m0,p_alp2,p_bet2);
	m0=mNaes1(length(ts1));
	[mNae,mNas,gamNa]=fm1(t,v,m0,p_alp2,p_bet2);

# Compute the estimated INa:
	m00=mNae.^3.*hNa.*(v-ENa);
	gNa=ci0mn/min(m00);
	INa=gNa*m00;
	INas1=gNa*mNaes1.^3.*hNas1.*(vs1-ENa);

	

endfunction
  
