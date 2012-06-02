/* C++ functions used in package "jmec" */

/* #include <math.h> */

#include <R.h>
#include <Rmath.h>

/* rough numeric estimation of mean and sd of candidate distribution for Q draws  */
/* extern "C" __declspec(dllexport) void cdest( */
extern "C" void cdest(
		double *lvT,
		double *lvC,
		int *nobsT,
		int *nobsC,
		int *nvarT,
		int *nvarC,
		double *cpT,
		double *cpC,
		double *betaTest,
		double *betaCest,
		double *thetaest,
		double *alphaest,
		double *bval,
		int *bvallength,
		double *lambda0Tv,
		double *lambda0Cv,
		double *lambda0Tvcum,
		double *lambda0Cvcum,
		int *eventv,
		int *censv,
		double *l,
		double *fb
		)
{
	int i, j;
	double const pi=4*atan(1);
	/* compute crossproduct linear predictor and betas for event */
		for(j = 1; j <= nobsT[0]; j++)
		{
			for (i = 1; i <= nvarT[0]; i++)
			{
				cpT[j - 1] = cpT[j - 1] + betaTest[i - 1]*lvT[j + (i - 1)*nobsT[0] - 1];
			};
		};
	/* compute crossproduct linear predictor and betas for censoring */
		for(j = 1; j <= nobsC[0]; j++)
		{
			for (i = 1; i <= nvarC[0]; i++)
			{
				cpC[j - 1] = cpC[j - 1] + betaCest[i - 1]*lvC[j + (i - 1)*nobsC[0] - 1];
			};
		};
	/* compute density for vector of b values */
		for(j = 0; j < bvallength[0]; j++)
		{
			for(i = 0; i < nobsT[0]; i++)
			{
				l[i] = 	pow(lambda0Tv[i]*exp(cpT[i] + bval[j]), eventv[i])*
						exp(-1*lambda0Tvcum[i]*exp(cpT[i] + bval[j]))*
						pow(lambda0Cv[i]*exp(cpC[i] + alphaest[0]*bval[j]), censv[i])*
						exp(-1*lambda0Cvcum[i]*exp(cpC[i] + alphaest[0]*bval[j]));
			};
			fb[j] = (1/sqrt(2*pi*thetaest[0]))*exp(-1*pow(bval[j],2)/(2*thetaest[0]))*l[0];
			for(i = 1; i < nobsT[0]; i++)
			{
				fb[j] = fb[j]*l[i];
			};
		};
}

/* Metropolis-Hastings sampler for Q draws  */
/* extern "C" __declspec(dllexport) void mh( */
extern "C" void mh(
		double *lvT,
		double *lvC,
		int *nobsT,
		int *nobsC,
		int *nvarT,
		int *nvarC,
		double *cpT,
		double *cpC,
		double *betaTest,
		double *betaCest,
		double *thetaest,
		double *alphaest,
		double *lambda0Tv,
		double *lambda0Cv,
		double *lambda0Tvcum,
		double *lambda0Cvcum,
		int *eventv,
		int *censv,
		double *l,
		double *bv,
		double *densv,
		double *bcanv,
		double *denscanv,
		double *av,
		double *rnorm1,
		double *walkv,
		double *runifv,
		int *Q
		)
{
	int i, j;
	double const pi=4*atan(1);
	/* compute crossproduct linear predictor and betas for event */
		for(j = 1; j <= nobsT[0]; j++)
		{
			for (i = 1; i <= nvarT[0]; i++)
			{
				cpT[j - 1] = cpT[j - 1] + betaTest[i - 1]*lvT[j + (i - 1)*nobsT[0] - 1];
			};
		};
	/* compute crossproduct linear predictor and betas for censoring */
		for(j = 1; j <= nobsC[0]; j++)
		{
			for (i = 1; i <= nvarC[0]; i++)
			{
				cpC[j - 1] = cpC[j - 1] + betaCest[i - 1]*lvC[j + (i - 1)*nobsC[0] - 1];
			};
		};
	/* Preparations for Metropolis-Hastings loop */
		bv[0] = rnorm1[0];
		for(i = 0; i < nobsT[0]; i++)
		{
			l[i] = 	pow(lambda0Tv[i]*exp(cpT[i] + bv[0]), eventv[i])*
					exp(-1*lambda0Tvcum[i]*exp(cpT[i] + bv[0]))*
					pow(lambda0Cv[i]*exp(cpC[i] + alphaest[0]*bv[0]), censv[i])*
					exp(-1*lambda0Cvcum[i]*exp(cpC[i] + alphaest[0]*bv[0]));
		};
		densv[0] = (1/sqrt(2*pi*thetaest[0]))*exp(-1*pow(bv[0],2)/(2*thetaest[0]))*l[0];
		for(i = 1; i < nobsT[0]; i++)
		{
			densv[0] = densv[0]*l[i];
		};
	/* Metropolis-Hastings loop */
		for(j = 1; j < Q[0]; j++)
		{
			bcanv[j] = bv[j-1] + walkv[j];
			for(i = 0; i < nobsT[0]; i++)
			{
				l[i] = 	pow(lambda0Tv[i]*exp(cpT[i] + bcanv[j]), eventv[i])*
						exp(-1*lambda0Tvcum[i]*exp(cpT[i] + bcanv[j]))*
						pow(lambda0Cv[i]*exp(cpC[i] + alphaest[0]*bcanv[j]), censv[i])*
						exp(-1*lambda0Cvcum[i]*exp(cpC[i] + alphaest[0]*bcanv[j]));
			};
			denscanv[j] = (1/sqrt(2*pi*thetaest[0]))*exp(-1*pow(bcanv[j],2)/(2*thetaest[0]))*l[0];
			for(i = 1; i < nobsT[0]; i++)
			{
				denscanv[j] = denscanv[j]*l[i];
			};
			av[j] = denscanv[j]/densv[j-1];
			if (av[j] >= 1)
			{
				bv[j] = bcanv[j];
				densv[j] = denscanv[j];
			}
			else
			{
				if(runifv[j] <= av[j])
				{
					bv[j] = bcanv[j];
					densv[j] = denscanv[j];
				}
				else
				{
					bv[j] = bv[j-1];
					densv[j] = densv[j-1];
				}
			}
		};
}
