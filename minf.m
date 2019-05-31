# THIS SOFTWARE COMES WITH NO WARRANTY WHATSOEVER, EXPRESSED OR IMPLIED. 
# USE IT AT YOUR OWN RISK!
#
#	By T.I. Toth, Cardiff University, U.K.; 1996-2002
#
#
#
# This is a simple Boltzmann curve
#
#	Input:
# v: voltage (can be a vector);
# q: "slope" par.;
# V0: "half value" par.;
#
#	Output:
# xm: value of the function at v;



function xm=minf(v,q,V0)

	xm=q*(v+V0);
	xm=1./(1+exp(xm));

endfunction	
