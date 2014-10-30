/*
 
 transfer_gadget.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 */

#include "transfer_function.hh"
#include <stdlib.h>
#include <string.h>

//! Implementation of abstract base class TransferFunction.
/*!
 This class implements the power spectrum files used in CMU simulations.
 the file lists log10 of k (in Mpc/h) and Pk( Mpc/h, ^-3).

 The tilt is included in the power spectrum; baryon and dm are the same.

 we convert back to transfer function from these files by removing k** ns.

 Yu Feng <yfeng1@berkeley.edu>

 Most of the core code is from MP-Gadget3/GenIC which originated from N-GenIC.

 */
class transfer_gadget_plugin : public transfer_function_plugin
{
protected:
	using transfer_function_plugin::cosmo_;
    std:: string m_filename_Pk;	
    struct pow_table {
        double logk;    
        double logd;    
    }* PowerTable;
    int NPowerTable;
    static int compare_logk(const void *a, const void *b) {
        if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk))
            return -1;

        if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk))
            return +1;

        return 0;
    }
public:
	//! Constructor for Eisenstein & Hu fitting for transfer function
	/*!
	 \param aCosm structure of type Cosmology carrying the cosmological parameters
	 \param Tcmb mean temperature of the CMB fluctuations (defaults to
	 Tcmb = 2.726 if not specified)
	 */
	transfer_gadget_plugin( config_file &cf )//Cosmology aCosm, double Tcmb = 2.726 )
    :  transfer_function_plugin(cf)
	{
		m_filename_Pk = pcf_->getValue<std::string>("cosmology","powerspectrum_file");
        read_power_table(m_filename_Pk.c_str());
	}

    void read_power_table(const char * FileWithInputSpectrum) {
        FILE *fd;
        char buf[500];
        double k, p;

        strcpy(buf, FileWithInputSpectrum);

        if(!(fd = fopen(buf, "r")))
        {
            printf("can't read input spectrum in file '%s'\n", buf);
            abort();
        }

        NPowerTable = 0;
        do
        {
            if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
                NPowerTable++;
            else
                break;
        }
        while(1);

        fclose(fd);

        printf("found %d pairs of values in input spectrum table\n", NPowerTable);
        fflush(stdout);

        PowerTable = (struct pow_table*)malloc(NPowerTable * sizeof(struct pow_table));

        int i = 0;
        sprintf(buf, FileWithInputSpectrum);

        if(!(fd = fopen(buf, "r")))
        {
            printf("can't read input spectrum in file '%s' \n", buf);
            abort();
        }

        i = 0;
        do
        {
            if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
            {
                PowerTable[i].logk = k;
                PowerTable[i].logd = p;
                i++;
            }
            else
                break;
        }
        while(1);

        fclose(fd);

        qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);
    }
	
    double PowerSpec_Tabulated(double k)
    {
        double logk, logd, P, kold, u, dlogk, Delta2;
        int binlow, binhigh, binmid;

        double mydlogk,dlogk_PowerTable;
        int mybinhigh,mybinlow,mybinmid;

        kold = k;

        logk = log10(k);

        if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
            return 0;

        dlogk_PowerTable = PowerTable[1].logk-PowerTable[0].logk;
        mydlogk = logk - PowerTable[0].logk;
        mybinlow = (int)(mydlogk/dlogk_PowerTable);
        mybinhigh = mybinlow+1;

        dlogk = PowerTable[mybinhigh].logk - PowerTable[mybinlow].logk;

        if(dlogk == 0)
            abort();

        u = (logk - PowerTable[mybinlow].logk) / dlogk;

        logd = (1 - u) * PowerTable[mybinlow].logd + u * PowerTable[mybinhigh].logd;

        P = pow(10.0, logd);
        return P;
    }
	//! Computes the transfer function for k in Mpc/h by calling TFfit_onek
	inline double compute( double k, tf_type type ){
        /* remove the tilt because this already has tilt! */
        return pow(PowerSpec_Tabulated(k) * pow((double)k, - (double)cosmo_.nspect), 0.5);
	}
	
	inline double get_kmin( void ){
		return 1e-4;
	}
	
	inline double get_kmax( void ){
		return 1.e4;
	}
	
};


namespace{
	transfer_function_plugin_creator_concrete< transfer_gadget_plugin > creator("gadget");
}

