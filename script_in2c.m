# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#	This is a driver script file in Octave for the function in_estim2c
#	and its outputs;
# N.B. The input data should be loaded or defined before using this programme.
#
#

1;		# This is to indicate that this is NOT a function

# Signal ready (beep):
	printf("\a\n");

# Edit the particulars (filename, k_range, name of the data matrix) below
# (file.name: name of the output file where the estimation results are saved;
#  k_range=[k1:k2], or k_range=[k1,k2,k3:k4,k5] etc.; actual name for the
#  data matrix burst should be used):

	fnam1="file.name"

	[fnum1,errmsg]=fopen(fnam1,"a");

# Print file header:
	fprintf(fnum1,"\n\n\n\n\n\n");
	fprintf(fnum1,"\t\t\t File: %s with parmH\n\n",fnam1);
	fprintf(fnum1,"    Estimtion interval (ms): [%4.1f,%5.1f]",tb,te);
	fprintf(fnum1,"    Stimulus start at %5.1f ms.\n\n\n",tsb);
	fprintf(fnum1," k0  Iin (pA)   Cm (pF)   gH (nS)   gL (nS)  EL (mV)  Er (mV)    Eh (mV)   Ih (pA)\n");

# Compute estimates using in_estim2c for every trace:
	for k0=k_range
	   [Cm,gH,gL,cih,Eh,EL,Er,ci2,t,v,vch,dv,IL,IH]=in_estim2c(burst,\
		tsb,tse,tb,te,k0,mch,dc_shift,cifc,parmH,EH);
	fprintf(fnum1,"%2d   %7.2f   %7.2f   %7.2f  %7.2f   %7.2f  %7.2f   %7.2f    %7.2f\n",k0,ci2,Cm,gH,gL,EL,Er,Eh,cih);
	   printf("rel_sq_err=%11.4e\n",norm(ci2-Cm*dv-IL-IH)/abs(ci2));
	endfor

fprintf(fnum1,"\n parmH: %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",parmH);

	fclose(fnum1);
	clear fnum1 fnam1
