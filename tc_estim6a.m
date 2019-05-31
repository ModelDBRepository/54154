# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This procedure finds an estimate for gK and for the parameters of
#	mK_inf(V), a simple Boltzmann curve, and gamma_K(V), a polynomial of
#	degree 4 in V.
#    
#    Notation: tb=t(1)=ts1(length(ts1)), te=t(length(t)), 
#	       tsb=ts1(1) : stim. start.
#
#
#		Input:
#   t: time vector on [tb,te];
#   v: V(t) on [tb,te];
#   ci0: ci0=IK+INa on [tb,te];
#   ts1: time from stim. start to AP;
#   vs1: corresponding V vector;
#   Er: resting pot.;
#   EK: K rev.pot.;
#   mKmx: assumed maximal value of mK(t) (in accordance with the
#	  steady-state activ. properties of IK); 0.79<mKmx<0.99 is used;
#   gNa: max. conductance of INa;
#   INa: INa(t) on [tb,te];
#   I0K: value of IK at Er, I0K<0.1 pA is chosen;
#
#
#		Output:
#   kb: index of the first `good' point of mK from which it 
#	is monotonic and positive up to its maximum; (corresponding to 
#	the time tkb);
#   ikmx: index where mK is maximal (has the largest positive peak);
#   ke: index of the last `good' point of mK at which its time 
#	derivative is still non-positive; (corresponding the time tke);
#   mKa: 'original' mK on [tb,te];
#   dmK: dmK/dt on [tkb,tke];
#   mch2: order+1 of Cheb. approx. used to compute dmK;
#   vg: the voltage vector on which gamma_K(V)>0; vg is a subset of v(kb:ke);
#   gamKa: gamma_K(V)) on vg;
#   gK: optimized value of gK=g0K^4;
#   qK: `slope' par. of mK_inf(V) (a Boltzmann curve);
#   V0K: `half value' of mK_inf(V);
#   mKes1: mK(t) on [tsb,tb];
#   mKss1: mK_inf(V(t)) on [tsb,tb];
#   gamKs1: gamma_K(V(t)) on [tsb,tb];
#   mKe: mK(t) on [tb,te];
#   mKs: mK_inf(V(t)) on [tb,te];
#   gamK: evaluated polynomial estimate of gamma_K(V(t)) on [tb,te];
#   p_gamK: vector of the polynom. coeffs. for gamma_K(V);
#   IKs1: reconstructed IK(t) on [tsb,tb];
#   IK: reconstructed IK(t) on [tb,te];
#   gNae: new, optimized value of gNa;
#   INae: new INa(t) on [tb,te], computed with gNae;
#   
#
#		Internal variables:
#   N1: length of the selected sub-interval [t(kb),t(ke)], N1=ke-kb+1;
#   a: a(t)=[IK/(v-EK)]^(1/4);
#   amx: amx=max(a)=a(ikmx);
#   vmx: vmx=v(ikmx);
#   mKEr: steady-state value of mK at Er;
#   g0K: g0K=gK^(1/4);
#   t1: t(kb:ke);
#   v1: v(kb:ke);
#   mch2x: maximal admissible order+1 of Cheb. approx.;
#   ng0: degree of the polynom. approx.;
#   G: coeff. matrix od the lsq. polyn. approx. for gammaK(V);
#   G1: coeff. matrix of the lin. lsq. proc. for gK and gNa;
#   sig, res: s.d. and residual vector of the lin. lsq.;


#		External procedures and functions:
#  fm2(): computes m(t), m_inf(V(t)), gamma(V(t)) when m_inf(V) is a
#	  Boltzmann curve, and gamma(V) is a polynomial in V.
#  minf(): Boltzmann curve;



function [kb,ikmx,ke,mKa,dmK,mch2,vg,gamKa,gK,qK,V0K,mKes1,mKss1,gamKs1,\
	mKe,mKs,gamK,p_gamK,IKs1,IK,gNae,INae]=\
	tc_estim6a(t,v,ci0,ts1,vs1,Er,EK,mKmx,gNa,INa,I0K)


	N0=length(t);		# no. of data points on [tb,te]
# Check consistency of the data:
	if ((length(v) !=N0)||(length(ci0)!=N0)||(length(INa)!=N0))
	   printf("data lengths: lt=%3d lv=%3d lci0=%3d lINa=%3d\n",\
		N0,length(v),length(ci0),length(INa));
	   error("In tc_estim6a: Data are inconsistent\n");
	endif

# Compute IK:
	IK=ci0-INa;

# Compute aK(t)=g0K*mK(t) and find its max.:
	a=max(IK./(v-EK),0);			# gK*mK^4
	a=a.^(1/4);
	amx=max(a);
	ikmx=max(find(1-sign(amx-a)));		# index of a_max
	for ii=ikmx-1:-1:1
	   if ((a(ii+1)<a(ii))||(a(ii+1)*a(ii)==0)) break; endif
	endfor
	kb=ii;			# start of the `usable' part of mK

# Find the last `usable' point of aK(t) after its max.:
	for ii=ikmx:N0-1
	   if (a(ii)<a(ii+1)) break; endif
	endfor
	ke=ii;
	t1=t(kb:ke);
	v1=v(kb:ke);
	N1=ke-kb+1;			# length of t1, v1 and dmK (below)
	ik1=ikmx-kb+1;
	ik2=ik1+1;
	ik0=ik1-1;

# Max. conductance:
	g0K=amx/mKmx;
	mKa=a/g0K;			# `original' mK data

# Estimate qK and V0K as parameters of a Boltzmann curve:
	vmx=v(ikmx);
	mKEr=(I0K/(Er-EK))^(1/4);
	mKEr=mKEr/g0K;
	qK=log((1/mKEr-1)/(1/mKmx-1))/(Er-vmx);
	V0K=log(1/mKmx-1)/qK-vmx;

# Compute dmK/dt over t1 such that its zero-crossing is at max mK:
	mch2x=10;
	for mch2=4:mch2x
	   dmK=df_ch_vect(t1,mKa(kb:ke),mch2)';
	   if ((dmK(ik1)*dmK(ik2)<=0) || (dmK(ik0)*dmK(ik1)<=0))
	      break;
	   endif
	endfor

# Compute gamma_K(V(t))=dmK/dt/(mK_inf(V(t))-mK(t)):
	gamKa=dmK./(minf(v1,qK,V0K)-mKa(kb:ke));
	gamKa(ik1)=0.5*(gamKa(ik0)+gamKa(ik2));		# average at mK_peak

# Exclude negative gamma_K values:
	kg=0;
	for k=1:N1
	   ga0=gamKa(k);
	   if (ga0>0)
	      kg++;
	      gamKa(kg)=ga0;
	      vg(kg)=v(kb+kg-1);
	   endif
	endfor
	gamKa=gamKa(1:kg);

# Estimate gamma_K(V) as a polynomial of degree 4 in V:
	ng0=4;
	clear G;
	G(:,ng0+1)=ones(size(vg));
	G(:,ng0)=vg;
	for ng=ng0:-1:2
	   G(:,ng-1)=G(:,ng).*vg;
	endfor
	[p_gamK,sig,res]=ols(gamKa,G);
	sig_gamK=sig
	max_err=max(abs(res));

# Compute mK(t), mK_inf(V(t)), and gamma_K(V(t)) during ts1 and t:	
	[mKes1,mKss1,gamKs1]=fm2(ts1,vs1,mKEr,qK,V0K,p_gamK);
	m0=mKes1(length(ts1));
	[mKe,mKs,gamK]=fm2(t,v,m0,qK,V0K,p_gamK);

# Optimize gK and gNa:
	G1(:,1)=mKe.^4.*(v-EK);
	G1(:,2)=INa/gNa;
	[gz,sig,res]=ols(ci0,G1);
	sig_curr=sig
	max_err_curr=max(abs(res))
	gK=gz(1);
	gNae=gz(2);
	INae=gNae*G1(:,2);
	IK=gK*G1(:,1);
	IKs1=gK*mKes1.^4.*(vs1-EK);


endfunction
  
