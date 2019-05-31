# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This procedure computes parameter estimates for INa, i.e. 
#	mNa_inf(V) as product of a steep Boltzmann curve and a polynomial
#	of degree ng1<9 in V, and gamma_Na(V) or 1/gamma_Na(V) as 
#	a polynomial of degree ng0<9 in V whichever yields the better fit.
#	Both gNa and gK are optimized by lin. lsq.
#	In this version (df), IK is increased such that an "a-loop" appears 
#	with no self-intersection in the V-a plane. This is ensured by 
#	checking that a(t) values on the descending phase of the AP 
#	are greater than the corresponding ones on the ascending phase.
#   Notation: tb=t(1)=ts1(length(ts1)), te=t(length(t)), 
#	      tsb=ts1(1): stim. start.
#
#
#		Input:
#   t: time vector on [tb,te];
#   v: V(t) on [tb,te];
#   ci0: ci0=IK+INa on [tb,te];
#   ts1: time from stim. start to AP: [tsb,tb];
#   vs1: corresponding V vector;
#   Er: resting pot.;
#   ENa: Na rev.pot.;
#   mNax0: (approx.) max. value of mNa_inf(V);
#   ng0: degree of the polynom. for gamma_Na(V) or 1/gamma_Na(V);
#   ng1: degree of the polynom. for maNa_inf(V);
#   hNa: hNa(t) on [tb,te];
#   hNas1: hNa(t) on [tsb,tb];
#   gK: max. conductance of IK;
#   IK: IK on [tb,te];
#
#
#		Output:
#   ivmx: index of the peak (max.) of v;
#   inmx: index of t=tmx where mNa is maximal (has the largest positive peak);
#   nne: index of the last point of the v-interval with v(nne) appr.=v(1);
#   mNaa: 'original' mNa on [tb,te];
#   vng: a subset of v(1:nne) where both mNa_inf(V) and gamma_Na(V) 
#	are computed and positive;
#   tng: time vector corresponding to vng;
#   mNasa: computed values of mNa_inf(V) on vng;
#   gamNaa: computed values of gamma_Na(V) on vng;
#   p_gamNa: coeff. vector of the polynom. estimate for gamma_Na(V)
#	or 1/gamma_Na(V); the last entry of p_gamNa shows which case
#	was accepted; +1: gamma_Na(V), -1: 1/gamma_Na(V);
#   p_mNa: coeff. vector of the polynom. estimate for mNa_inf(V);
#   q1: `slope' factor of the Boltzmann curve of mNa_inf(V);
#   V0: Vhalf of the Boltzmann curve of mNa_inf(V);
#   gNa: optimized value of gNa;
#   da: da/dt where a=g0Na*mNa(t);
#   mch2: order+1 of the Cheb. approx used to compute da;
#   mNaes1: mNa(t) computed on [tsb,tb];
#   mNass1: mNa_inf(V(t)) computed on [tsb,tb];
#   gamNas1: gamma_Na(V(t)) computed on [tsb,tb];
#   mNae: mNa(t) computed on [tb,te];
#   mNas: mNa_inf(V(t)) computed on [tb,te];
#   gamNa: gamma_Na(V(t)) computed on [tb,te];
#   INas1: reconstructed INa(t) on [tsb,tb];
#   INa: reconstructed INa(t) on [tb,te];
#   gKe: new, optimized value of gK;
#   IKe: new IK on [tb,te] computed with gKe;
#   norm2_curr: quadratic error of ci0-IKe-INa;
#   final_sig: min. sigma of the fit of gNa and gKe;
#   final_err: maximal error of the fit; 
#   
#
#		Internal variables:
#   ci0mx1: limiting current value for IK;
#   tvmx: t(ivmx) where V is maximal;
#   bK: mK^4*(v-EK);
#   dt1: sampling time interval, dt1=t(2)-t(1);
#   thst: (varied) min. distance between t(inmx) and t(ivmx);
#   thstmx: max. value of thst;
#   cf: factor used to increase IK;
#   dcf: (preset) increment of cf;
#   ramn: relative min. distance of corresponding points in the "a-loop";
#   atol: threshold for ramn;
#   sig1: preceding error value for comparison with the actual one; 
#   g0Na: g0Na=gNa^(1/3);
#   a: g0Na*mNa;
#   mch2x: max. admissible order+1 of Cheb. approx.;
#   t1: t(1:ivmx-1);
#   v1: v(1:ivmx-1);
#   t2: t(nne:-1:ivmx+1);
#   v2: v(nne:-1:ivmx+1);
#   a1,a2,da1,da2: sub-vectors of a and da corresponding to t1 and t2, resp.;
#   Ha: coeff. matrix of the lin. eq.sys. for gamma(V) and g0Na*alpha(V)
#	at each pair of identical V-values from v1 and v2;
#   aNa0: g0Na*mNa_inf(V) as obtained from the solution to the lin. eq.sys.;
#   G: coeff. matrix of the lsq. polyn. approx. for mNa_inf(V) and gamma_Na(V);
#   G1: coeff. matrix of the lin. least square proc. for gNa and gK;
#   sig, res: sigma and residual vector of the least square proc.s;
#   inmx0,vng0,tng0 etc. : store the last `best' results;
#
#
#		External procedures and functions:
#  fm2df(): computes m(t), m_inf(V(t)), and gamma(V(t)) when m_inf(V) is
#	   a product of a Boltzmann curve and a polynomial in V, and
#	   gamma(V) or 1/gamma(V) is a polynomial in V.
#  minf(): Boltzmnann curve;




function [ivmx,inmx,nne,mNaa,vng,tng,mNasa,gamNaa,p_gamNa,p_mNa,q1,V0,gNa,\
	da,mch2,mNaes1,mNass1,gamNas1,mNae,mNas,gamNa,INas1,INa,gKe,\
	IKe,norm2_curr,final_sig,final_err]=\
	tc_estim7df(t,v,ci0,ts1,vs1,Er,ENa,mNax0,ng0,ng1,hNa,hNas1,gK,IK)


	N0=length(t);		# no. of data points on [tb,te]
# Check consistency of the data:
	if ((length(v) !=N0)||(length(ci0)!=N0)||(length(IK)!=N0))
	   printf("data lengths: lt=%3d lv=%3d lci0=%3d  lIK=%3d\n",\
		N0,length(v),length(ci0),length(IK));
	   error("In tc_estim7df: Data are inconsistent\n");
	endif

# Set and compute some test values for comparisons:
	ci0mx1=2.0*max(abs(ci0)); 			# bound for IK
	ivmx=max(find(1-sign(max(v)-v)));		# index of Vmax
	tvmx=t(ivmx);					# where V is maximal
	inmx=ivmx;
	inmx0=0;			# no useful results yet
	bK=IK/gK;			# mK^4*(v-EK)
	atol=0.02;			# minimal relative distance
	dcf=0.02;			# increment of cf (see below);

# Split the fct. V(t) into two monotonic (ascend. and descend.) parts:
	V1=v(1);
	nne=min(find(1+sign(V1-v(ivmx:N0))))+ivmx-1;
	v1=v(1:ivmx-1);
	t1=t(1:ivmx-1);
	v2=v(nne:-1:ivmx+1);
	t2=t(nne:-1:ivmx+1);

	dt1=t(2)-t(1);
	thst=2*dt1;		# separation time between t_vmax and t_mNa_max
	thst0=0;
	thstmx=6*dt1+1e-6;
	sig1=1e20;
	final_sig=sig1;
	final_err=sig1;
# Start a `big loop' which chooses the optimal thst=t(inmx)-t(ivmx)>0
    while (thst<thstmx)			# start of the `thst' loop
	cf=1;				# initial values
	ramn=0;
# Test whether the peaks of V and mNa are well separated (Vmax must come first)
# and whether there is a proper "a-loop" (see above):
	while ((t(inmx)-tvmx<thst)||(ramn<atol))
# Compute INa:
	   INa=ci0-cf*IK;
	   a=max(INa./(hNa.*(v-ENa)),0);
	   a=a.^(1/3);
# Find the largest positive peak of a:
	   mxa=max(a);
	   inmx=find(1-sign(mxa-a));		# index of mNa_max
# Exit loop if IK grows too large:
	   if (ci0mx1<max(cf*IK)) break;  endif
	   cf=cf+dcf;				# increase factor of IK
	   a1=a(1:ivmx-1);
	   a2=a(nne:-1:ivmx+1);
	   for ii=1:ivmx-1
	      vb=v1(ii);
	      a1a=a1(ii);
	      jj=0;
	      jj=min(find(1+sign(v2-vb)));
	      if (jj>1)
	         v2b=v2(jj);
	         a2b=a2(jj);
	         jj--;
	         v2a=v2(jj);
	         a2a=a2(jj);
	      else
		 ra(ii)=1e19;
	         continue;
	      endif
	      a1b=(vb-v2a)*(a2b-a2a)/(v2b-v2a)+a2a;
	      ra(ii)=(a1b-a1a)/mxa;
	   endfor
	   ramn=min(ra);
	endwhile		# end of ramn<atol loop

	if (ci0mx1<max(cf*IK))
	   if (inmx0==0)			# no useful results yet
	      printf("max V at %5.2f  max mNa at %5.2f\n",t(ivmx),t(inmx));
	      printf("Inactivation of INa is conflicting with the data.\n");
	      printf("max_ci0=%13.5g  max_IK=%13.5g  min_INa=%13.5g\n",
		max(ci0),max(IK),min(INa));
	      printf("In tc_estim7df: IK is too large\n");
	   endif
	   break;				# exit `thst' loop, too!
	endif

# Compute the actual diff. between the peaks of mNa and V:
	thst=t(inmx)-tvmx;

# Find the corresponding interpolated values of a(t(V)) and da(t(V)):
	mch2x=10;
	for mch2=4:mch2x
	   da=df_ch_vect(t,a,mch2)';		# da/dt
	   if (da(inmx-1)*da(inmx+1)<0) break; endif
	endfor
	da1=da(1:ivmx-1);			# da/dt on t1
	da2=da(nne:-1:ivmx+1);			# da/dt on t2

	ng=0;
	for ii=1:ivmx-1
	   vb=v1(ii);
	   a1a=a1(ii);
	   da1a=da1(ii);
	   jj=0;
	   jj=min(find(1+sign(v2-vb)));
	   if (jj>1)
	      v2b=v2(jj);
	      a2b=a2(jj);
	      da2b=da2(jj);
	      jj--;
	      v2a=v2(jj);
	      a2a=a2(jj);
	      da2a=da2(jj);
	   else
	      ii
	      continue;
	   endif
	   a1b=(vb-v2a)*(a2b-a2a)/(v2b-v2a)+a2a;
	   da1b=(vb-v2a)*(da2b-da2a)/(v2b-v2a)+da2a;
	   Ha=[1,-a1a;1,-a1b];
	   a_gam=Ha\[da1a;da1b];
	   gam=a_gam(2);
# Allow only positive solutions:
	   if ((gam>0)&&(a_gam(1)>0))
	      ng++;
	      vng(ng)=vb;
	      tng(ng)=t(ii); 
	      aNa0(ng)=a_gam(1)/gam;
	      gamNaa(ng)=gam;
	   endif
	endfor

	for ii=ivmx+1:nne
	   vb=v(ii);
	   a2a=a(ii);
	   da2a=da(ii);
	   jj=0;
	   jj=min(find(1+sign(v1-vb)));
	   if (jj>1)
	      v1b=v1(jj);
	      a1b=a1(jj);
	      da1b=da1(jj);
	      jj--;
	      v1a=v1(jj);
	      a1a=a1(jj);
	      da1a=da1(jj);
	   else
	      ii
	      continue;
	   endif
	   a2b=(vb-v1a)*(a1b-a1a)/(v1b-v1a)+a1a;
	   da2b=(vb-v1a)*(da1b-da1a)/(v1b-v1a)+da1a;
	   Ha=[1,-a2a;1,-a2b];
	   a_gam=Ha\[da2a;da2b];
	   gam=a_gam(2);
# Allow only positive solutions:
	   if ((gam>0)&&(a_gam(1)>0))
	      ng++;
	      vng(ng)=vb;
	      tng(ng)=t(ii);
	      aNa0(ng)=a_gam(1)/gam;
	      gamNaa(ng)=gam;
	   endif
	endfor

# Estimate gamma(V) or 1/gamma_Na(V) as a polynomial of degree ng0<9 in V:
# Usually ng0=8 or ng0=7;
	clear G;
	G(:,ng0+1)=ones(size(vng));
	G(:,ng0)=vng;
	for ng=ng0:-1:2
	   G(:,ng-1)=G(:,ng).*vng;
	endfor
# gamma(V):
	[p_gamNa1,sig,res]=ols(gamNaa,G);
	 p_gamNa1=[p_gamNa1;1];
# 1/gamma(V):
	[p_gamNa2,sig,res]=ols(1./gamNaa,G);
	 p_gamNa2=[p_gamNa2;-1];

# Compute g0Na and mNa:
	g0Na=max(aNa0)/mNax0;	# max(mNa_inf(V)) is approx. mNax0
	mNaa=a/g0Na;		# `original' mNa(t)
	mNasa=aNa0/g0Na;	# computed mNa_inf(V) on vng

# Estimate par.s of the steep Boltzmann curve of mNa_inf(V):
	V1=v(1);
	V2=v(4);
	q1=log((1/mNax0-1)/(1/mNasa(1)-1))/(V2-V1);
	V0=log(1/mNasa(1)-1)/q1-V1;
# Estimate mNa_inf(V) as a polynomial of degree <9 in V:
# Usually ng1=8 or ng1=7;
	clear G;
	G(:,ng1+1)=ones(size(vng));
	G(:,ng1)=vng;
	for ng=ng1:-1:2
	   G(:,ng-1)=G(:,ng).*vng;
	endfor
	[p_mNa,sig,res]=ols(mNasa,G);
	sig_mNa=sig
	max_err=max(abs(res));

# Compute mNa(t), mNa_inf(V(t)), and gamma_Na(V(t)) on [tsb,tb] and [tb,te]
# with p_gamNa1:
	bad_idx1=0;
	[mNaes1,mNass1,gamNas1]=fm2df(ts1,vs1,0.,p_mNa,q1,V0,p_gamNa1);
	m0=mNaes1(length(ts1));
	[mNae,mNas,gamNa]=fm2df(t,v,m0,p_mNa,q1,V0,p_gamNa1);
# Estimate gNa and gK if possible:
	if ((max(mNae)>=1.0)||(min(mNae)<0.)||\
		(max(mNas)>=1.0)||(min(mNas)<0.)||(min(gamNa)<0.))
	   bad_idx1=1;
	   siga=1e18;
	else
# Optimize gK and gNa:
	   clear G1;
	   G1(:,1)=bK;
	   G1(:,2)=mNae.^3.*hNa.*(v-ENa);
	   [gz1,siga,res]=ols(ci0,G1);
	   siga
	   max_erra=max(abs(res))
	endif

# Compute mNa(t), mNa_inf(V(t)), and gamma_Na(V(t)) on [tsb,tb] and [tb,te]
# with p_gamNa2:
	bad_idx2=0;
	[mNaes1,mNass1,gamNas1]=fm2df(ts1,vs1,0.,p_mNa,q1,V0,p_gamNa2);
	m0=mNaes1(length(ts1));
	[mNae,mNas,gamNa]=fm2df(t,v,m0,p_mNa,q1,V0,p_gamNa2);
# Estimate gNa and gK if possible:
	if ((max(mNae)>=1.0)||(min(mNae)<0.)||\
		(max(mNas)>=1.0)||(min(mNas)<0.)||(min(gamNa)<0.))
	   bad_idx2=1;
	   sigb=1e18;
	else
# Optimize gK and gNa:
	   clear G1;
	   G1(:,1)=bK;
	   G1(:,2)=mNae.^3.*hNa.*(v-ENa);
	   [gz2,sigb,res]=ols(ci0,G1);
	   sigb
	   max_errb=max(abs(res))
	endif

# Ignore meaningless results or save good ones:
	if ((bad_idx1)&&(bad_idx2))
	   thst=thst+dt1;
	   continue;
	elseif (siga<sigb)
	    sig=siga;
	    max_err=max_erra;
	    gz=gz1;
	    p_gamNa=p_gamNa1;
	else
	    sig=sigb;
	    max_err=max_errb;
	    gz=gz2;
	    p_gamNa=p_gamNa2;
	endif

	if (sig1>sig)
	   inmx0=inmx;
	   vng0=vng;
	   tng0=tng;
	   gamNaa0=gamNaa;
	   p_gamNa0=p_gamNa;
	   mNaa0=mNaa;
	   mNasa0=mNasa;
	   p_mNa0=p_mNa;
	   q10=q1;
	   V00=V0;
	   gz0=gz;
	   final_err=max_err;
	   final_sig=sig;
	   thst0=thst;
	   da0=da;
	   mch20=mch2;
	   sig1=sig;
	endif
	thst=thst+dt1;
    endwhile 			# end of the `thst' loop

	thst0
	final_sig
	final_err
# Set outputs to zero if no meaningful result:
	if (thst0==0) 
	   mNaa=0;vng=0;tng=0;mNasa=0;gamNaa=0;p_gamNa=0;p_mNa=0;
	   q1=0;V0=0;da=0;mch2=0;mNaes1=0;mNass1=0;gamNas1=0;
	   mNae=0;mNas=0;gamNa=0;
	   norm2_curr=final_err;
	   gKe=0; gNa=0; INa=0; IKe=0; INas1=0;
	   printf("No meaningful results!\n");
	   return;
	endif
	
# Restore the `optimal' values:
	gKe=gz0(1);
	gNa=gz0(2);
	inmx=inmx0;
	vng=vng0;
	tng=tng0;
	gamNaa=gamNaa0;
	p_gamNa=p_gamNa0;
	mNaa=mNaa0;
	mNasa=mNasa0;
	p_mNa=p_mNa0;
	q1=q10;
	V0=V00;
	da=da0;
	mch2=mch20;

# Compute final mNa(t), mNa_inf(V(t)),gamma_Na(V(t)), INa(t) on [tsb,tb] 
# and [tb,te] with best p_gamNa:
	[mNaes1,mNass1,gamNas1]=fm2df(ts1,vs1,0.,p_mNa,q1,V0,p_gamNa);
	m0=mNaes1(length(ts1));
	[mNae,mNas,gamNa]=fm2df(t,v,m0,p_mNa,q1,V0,p_gamNa);
	INa=gNa*mNae.^3.*hNa.*(v-ENa);
	IKe=gKe*bK;
	INas1=gNa*mNaes1.^3.*hNas1.*(vs1-ENa);
# Quadratic error:
	norm2_curr=norm(ci0-IKe-INa)

	
endfunction
  
