# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
#	This is a driver script file in Octave for the functions tc_estim6a,
#	and tc_estim7df to find the `best value' of mKmx, the max of mK(t);


1;		# This is to indicate that this is NOT a function


# Signal ready:
	printf("\a\n");

	clear;
# Initial value and increment and upper limit of mKmx:
	mKmx=0.8;
	dm=0.015;
	mKmx0=0.986;


# Edit the particulars (filename etc) below:

	fnam1="new_estimation.res"	
	fver0="6a_7df"

	[fnum1,errmsg]=fopen(fnam1,"w");

# Print file header:
	fprintf(fnum1,"\t\t File: %s   version: %s\n\n",fnam1,fver0);
fprintf(fnum1,"\t mKmx \t      sigma/10^6      max_err/10^3     norm2_min/10^3\n");


# Run tc_estim6a and tc_estim7df with a series of mKmx values :
# Edit the name of the file to be loaded as required
	while (mKmx<mKmx0)
	   load -force input_data.dat
	   mKmx
	   [kb,ikmx,ke,mKa,dmK,mch2,vg,gamKa,gK,qK,V0K,mKes1,mKss1,gamKs1,\
		mKe,mKs,gamK,p_gamK,IKs1,IK,gNae,INae]=\
		tc_estim6a(t,v,ci0,ts1,vs1,Er,EK,mKmx,gNa,INa,I0K);
	   norm2_curr=norm(ci0-IK-INae)

	   ng0=8;
	   ng1=8;
	     [ivmx,inmx,nne,mNaa,vng,tng,mNasa,gamNaa,p_gamNa,p_mNa,q1,V0,gNa,\
		da,mch2,mNaes1,mNass1,gamNas1,mNae,mNas,gamNa,INas1,INa,gKe,\
		IKe,norm2_curr,final_sig,final_err]=tc_estim7df(t,v,ci0,\
		ts1,vs1,Er,ENa,mNax0,ng0,ng1,hNa,hNas1,gK,IK);
	   fprintf(fnum1,"%13.6f   %13.6g    %13.6g     %13.6g \n",\
		mKmx,1.e-6*final_sig,1.e-3*final_err,1.e-3*norm2_curr);
	   mKmx=mKmx+dm;
	endwhile

	fclose(fnum1);

	clear fnum1 fnam1 fver0
