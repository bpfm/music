/*
 
 output_art.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2012  Jose Onorbe & Oliver Hahn
 
 */
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include "output.hh"

template<typename T>
inline T bytereorder(T v )
{
	T rval;
	(reinterpret_cast<unsigned char*>(&rval))[3] = (reinterpret_cast<unsigned char*>(&v))[0];
	(reinterpret_cast<unsigned char*>(&rval))[2] = (reinterpret_cast<unsigned char*>(&v))[1];
	(reinterpret_cast<unsigned char*>(&rval))[1] = (reinterpret_cast<unsigned char*>(&v))[2];
	(reinterpret_cast<unsigned char*>(&rval))[0] = (reinterpret_cast<unsigned char*>(&v))[3];
	return rval;
}


template< typename T_store=float >
class art_output_plugin : public output_plugin
{
public:
	bool do_baryons_;
    bool swap_endianness_;
	double omegab_, omegam_;
	double gamma_;
    double astart_;
    double zstart_;
	size_t npcdm_;
    int hsize_;
    
protected:
    
    enum iofields {
		id_dm_mass, id_dm_vel, id_dm_pos
    };

	typedef struct io_header
	{
		char head[45];
		float aexpN; // current expansion factor
       	float aexp0; // initial expansion factor
		float amplt; // Amplitude of density fluctuations
		float astep; // Delta a -> time step. 
				// This value is also stored in pt.dat (binary 1 float)
				// It is recalculated by art for the next steps so just a small value should work
		int istep; // step (=0 in IC)
		float partw; // mass of highest res particle.
        float TINTG; //=0 in IC
        float EKIN; //SUM 0.5 * m_i*(v_i**2) in code units
        float EKIN1; //=0 in IC
        float EKIN2; //=0 in IC
        float AU0; //=0 in IC
        float AEU0; //=0 in IC
        int NROWC; // Number of particles in 1 dim (number of particles per page = NROW**2) 
	    int NGRIDC; // Number of cells in 1 dim
        int nspecies; // number of dm species
	    int Nseed; // random number used ( 0 for MUSIC? or set the random number used in the lowest level?)
        float Om0; //Omega_m
	    float Oml0; //Omega_L
        float hubble; //hubble constant h=H/100
	    float Wp5; // 
        float Ocurv; //Omega_k
	    float Omb0; // this parameter only appears in header in hydro runs
		float wpart[10]; // extras[0-9] particle masses from high res to low res (normalized to low res particle)
		    //  Mass of smallest particle=wpart[0]*0m0*2.746e+11*(Box/NGRID)**3 -> Msun/h
		    //  Mass of largest  particle=wpart[nspecies-1]*0m0*2.746e+11*(Box/NGRID)**3 -> Msun/h
		int lpart[10]; // extras[10-19] number of particles from high res to low res cumulative!!! 
		    //(i.e., lpart[0]=Nhigh res particles; lpart[1]=lpart[0]+N_this_level; etc) so lpart[nspecies-1]=N total
	    float extras[80]; //extras[20-99] 
		     //extras[9]=iLblock ->0 in IC 
             //extras[10]=LevMin  ->0 in IC
             //extras[11]=LevSmall ->0 in IC
             //extras[12]=LevLarge ->0 in IC
             //extras[13]=Omegab  ->0 in IC; fix it?
             //extras[14]=sig8    ->0 in IC; fix it?
             //extras[15]=Spslope ->0 in IC; fix it? Slope of the Power spectrum
             //extras[16]=iDEswtch ->0 in IC; DE Flag=0:LCDM 1:w 2:RP 3:SUGRA
             //extras[17]=DEw0    ->0 in IC; w0 for DE z=0
             //extras[18]=DEwprime ->0 in IC; DE parameter
		     //extras[59]= 0 or 1; is used as switch for random numbers generators [do not apply in music use 0?]
		     //extras[60]= lux - level of luxury  [do not apply in music use 0?]
		     //extras[79]=Lbox (Mpc/h)

	}header;

	typedef struct io_ptf
	{
		float astep;
	}ptf;
	
	header header_;
	ptf ptf_;
	std::string fname;
	size_t np_fine_gas_, np_fine_dm_, np_coarse_dm_;
	size_t block_buf_size_;
	size_t npartmax_;
    
    double YHe_;
	
	
	// helper class to read temp files
	class pistream : public std::ifstream
	{
	public:
		pistream (std::string fname, size_t npart )
		: std::ifstream( fname.c_str(), std::ios::binary )
		{
			size_t blk;
			
			if( !this->good() )
			{	
				LOGERR("Could not open buffer file in ART output plug-in");
				throw std::runtime_error("Could not open buffer file in ART output plug-in");
			}
			
			this->read( (char*)&blk, sizeof(size_t) );
			
			if( blk != (size_t)(npart*sizeof(T_store)) )
			{	
				LOGERR("Internal consistency error in ART output plug-in");
				LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
				throw std::runtime_error("Internal consistency error in ART output plug-in");
			}
		}
		
		pistream ()
		{
			
		}
		
		void open(std::string fname, size_t npart )
		{
			std::ifstream::open( fname.c_str(), std::ios::binary );
			size_t blk;
			
			if( !this->good() )
			{	
				LOGERR("Could not open buffer file \'%s\' in ART output plug-in",fname.c_str());
				throw std::runtime_error("Could not open buffer file in ART output plug-in");
			}
			
			this->read( (char*)&blk, sizeof(size_t) );
			
			if( blk != (size_t)(npart*sizeof(T_store)) )
			{	
				LOGERR("Internal consistency error in ART output plug-in");
				LOGERR("Expected %d bytes in temp file but found %d",npart*(unsigned)sizeof(T_store),blk);
				throw std::runtime_error("Internal consistency error in ART output plug-in");
			}
		}
	};
	
	
	// non-public member functions
	void write_header_file( void ) //PMcrd.DAT
	{
        std::string partfname = fname_ + "/PMcrd.dat";
        std::ofstream ofs( partfname.c_str(), std::ios::trunc );
        //ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
		header this_header(header_);
		//int blksize = sizeof(header); 
        //Should be 529 in a dm only run; 533 in a baryon run
        //but not working for alignment so:
        int blksize = hsize_;
		ofs.write( (char *)&blksize, sizeof(int) );
		//ofs.write( (char *)&this_header,sizeof(header));  //Not working because struct aligment, so:
		ofs.write( (char *)&this_header.head,sizeof(this_header.head));   
		ofs.write( (char *)&this_header.aexpN,sizeof(this_header.aexpN));   
		ofs.write( (char *)&this_header.aexp0,sizeof(this_header.aexp0));   
		ofs.write( (char *)&this_header.amplt,sizeof(this_header.amplt));   
		ofs.write( (char *)&this_header.astep,sizeof(this_header.astep));   
		ofs.write( (char *)&this_header.istep,sizeof(this_header.istep));   
		ofs.write( (char *)&this_header.partw,sizeof(this_header.partw));   
		ofs.write( (char *)&this_header.TINTG,sizeof(this_header.TINTG));   
		ofs.write( (char *)&this_header.EKIN,sizeof(this_header.EKIN));   
		ofs.write( (char *)&this_header.EKIN1,sizeof(this_header.EKIN1));   
		ofs.write( (char *)&this_header.EKIN2,sizeof(this_header.EKIN2));   
		ofs.write( (char *)&this_header.AU0,sizeof(this_header.AU0));   
		ofs.write( (char *)&this_header.AEU0,sizeof(this_header.AEU0));   
		ofs.write( (char *)&this_header.NROWC,sizeof(this_header.NROWC));   
		ofs.write( (char *)&this_header.NGRIDC,sizeof(this_header.NGRIDC));   
		ofs.write( (char *)&this_header.nspecies,sizeof(this_header.nspecies));   
		ofs.write( (char *)&this_header.Nseed,sizeof(this_header.Nseed));   
		ofs.write( (char *)&this_header.Om0,sizeof(this_header.Om0));   
		ofs.write( (char *)&this_header.Oml0,sizeof(this_header.Oml0));   
		ofs.write( (char *)&this_header.hubble,sizeof(this_header.hubble));   
		ofs.write( (char *)&this_header.Wp5,sizeof(this_header.Wp5));   
		ofs.write( (char *)&this_header.Ocurv,sizeof(this_header.Ocurv));   
		if (do_baryons_) {
		    ofs.write( (char *)&this_header.Omb0,sizeof(this_header.Omb0));   
        }
		ofs.write( (char *)&this_header.wpart,sizeof(this_header.wpart));   
		ofs.write( (char *)&this_header.lpart,sizeof(this_header.lpart));   
		ofs.write( (char *)&this_header.extras,sizeof(this_header.extras));   
		ofs.write( (char *)&blksize, sizeof(int) );
		ofs.close();
		LOGINFO("ART : done writing header file.");
	}

	void write_pt_file( void ) //pt.dat
	{
        std::string partfname = fname_ + "/pt.dat";
        std::ofstream ofs( partfname.c_str(), std::ios::trunc );
        //ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );
		ptf this_ptf(ptf_);
		int blksize = sizeof(ptf); //4
		ofs.write( (char *)&blksize, sizeof(int) );
		ofs.write( (char *)&this_ptf,sizeof(ptf));
		ofs.write( (char *)&blksize, sizeof(int) );
		ofs.close();
		LOGINFO("ART : done writing pt file.");
	}
    
    
    void adjust_buf_endianness( T_store* buf )
    {
        if( swap_endianness_ )
        {
            for( size_t i=0; i<block_buf_size_; ++i )
                buf[i] = bytereorder<T_store>( buf[i] );
        }
    }

	/*
     The direct format write the particle data in pages. Each page of particles is read into a common block,
	 which has the structure: X(Npage),Y(Npage),Z(Npage),Vx(Npage),Vy(Npage),Vz(Npage). 
	 The number of particles in each page (Npage) is Npage = Nrow**2; Npages = (N_particles -1)/NPAGE +1
	 so in last page sometimes can be tricky: N_in_last=N_particles -NPAGE*(Npages-1)

	 There are NO Fortran size blocks pre or after these blocks!!

	 Contradiction with documentation?? one file for each type of particle
     however Daniel sent me just one file for a zoom all particle info together. 
    */
	void assemble_DM_file( void ) //PMcrs0.dat
	{
		// have to fix file name
		std::string partfname = fname_ + "/PMcrs0.dat";
		std::ofstream ofs( partfname.c_str(), std::ios::trunc );
		
		// generate all temp file names
		char fnx[256],fny[256],fnz[256],fnvx[256],fnvy[256],fnvz[256],fnm[256];
		sprintf( fnx,  "___ic_temp_%05d.bin", 100*id_dm_pos+0 );
		sprintf( fny,  "___ic_temp_%05d.bin", 100*id_dm_pos+1 );
		sprintf( fnz,  "___ic_temp_%05d.bin", 100*id_dm_pos+2 );
		sprintf( fnvx, "___ic_temp_%05d.bin", 100*id_dm_vel+0 );
		sprintf( fnvy, "___ic_temp_%05d.bin", 100*id_dm_vel+1 );
		sprintf( fnvz, "___ic_temp_%05d.bin", 100*id_dm_vel+2 );
		sprintf( fnm,  "___ic_temp_%05d.bin", 100*id_dm_mass  );
		
		// create buffers for temporary data
		T_store *tmp1, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6, *tmp7;
		
		tmp1 = new T_store[block_buf_size_];
		tmp2 = new T_store[block_buf_size_];
		tmp3 = new T_store[block_buf_size_];
		tmp4 = new T_store[block_buf_size_];
		tmp5 = new T_store[block_buf_size_];
		tmp6 = new T_store[block_buf_size_];
		tmp7 = new T_store[block_buf_size_];
		
		
		// read in the data from the temporary files in slabs and write it to the output file
		size_t npleft, n2read;
		size_t npcdm = npcdm_;
		
		LOGINFO("writing DM data to ART format file");
        //ofs.open(fname_.c_str(), std::ios::binary|std::ios::trunc );

		pistream ifs_x, ifs_y, ifs_z, ifs_vx, ifs_vy, ifs_vz, ifs_m;
		
		ifs_x.open( fnx, npcdm );
		ifs_y.open( fny, npcdm );
		ifs_z.open( fnz, npcdm );
		ifs_vx.open( fnvx, npcdm );
		ifs_vy.open( fnvy, npcdm );
		ifs_vz.open( fnvz, npcdm );
		ifs_m.open( fnm, npcdm );
		
		npleft = npcdm;
		n2read = std::min(block_buf_size_,npleft);
		while( n2read > 0 )
		{
			ifs_x.read( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
			ifs_y.read( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
			ifs_z.read( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
			ifs_vx.read( reinterpret_cast<char*>(&tmp4[0]), n2read*sizeof(T_store) );
			ifs_vy.read( reinterpret_cast<char*>(&tmp5[0]), n2read*sizeof(T_store) );
			ifs_vz.read( reinterpret_cast<char*>(&tmp6[0]), n2read*sizeof(T_store) );
			ifs_m.read( reinterpret_cast<char*>(&tmp7[0]), n2read*sizeof(T_store) );
            
            adjust_buf_endianness( tmp1 );
            adjust_buf_endianness( tmp2 );
            adjust_buf_endianness( tmp3 );
            adjust_buf_endianness( tmp4 );
            adjust_buf_endianness( tmp5 );
            adjust_buf_endianness( tmp6 );
            adjust_buf_endianness( tmp7 );
		
            ofs.write( reinterpret_cast<char*>(&tmp1[0]), n2read*sizeof(T_store) );
            ofs.write( reinterpret_cast<char*>(&tmp2[0]), n2read*sizeof(T_store) );
            ofs.write( reinterpret_cast<char*>(&tmp3[0]), n2read*sizeof(T_store) );
            ofs.write( reinterpret_cast<char*>(&tmp4[0]), n2read*sizeof(T_store) );
            ofs.write( reinterpret_cast<char*>(&tmp5[0]), n2read*sizeof(T_store) );
            ofs.write( reinterpret_cast<char*>(&tmp6[0]), n2read*sizeof(T_store) );

			npleft -= n2read;
			n2read = std::min( block_buf_size_,npleft );
		}
		
		ifs_x.close();
		ifs_y.close();
		ifs_z.close();
		ifs_vx.close();
		ifs_vy.close();
		ifs_vz.close();
		ifs_m.close();
		ofs.close();
		
		// clean up temp files
        unlink(fnx);
		unlink(fny);
		unlink(fnz);
		unlink(fnvx);
		unlink(fnvy);
		unlink(fnvz);
		unlink(fnm);

        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] tmp4;
        delete[] tmp5;
        delete[] tmp6;
        delete[] tmp7;
		
		LOGINFO("ART : done writing DM file.");
		
	}
	
	

public:


	explicit art_output_plugin ( config_file& cf )
	: output_plugin( cf )
	{
		if( mkdir( fname_.c_str(), 0777 ) )
                {
                        perror( fname_.c_str() );
                        throw std::runtime_error("Error in art_output_plugin!");
                }

		do_baryons_ = cf.getValueSafe<bool>("setup","baryons",false);
        // fix to header size (alignment problem)
		if (!do_baryons_)
            hsize_ = 529; // dm run
		else
            hsize_ = 533; // hydro run

		omegab_  = cf.getValueSafe<double>("cosmology","Omega_b",0.0);
	    omegam_  = cf.getValue<double>("cosmology","Omega_m");
        zstart_  = cf.getValue<double>("setup","zstart");
        astart_ = 1.0/(1.0+zstart_);
        
        
        swap_endianness_ = cf.getValueSafe<bool>("output","art_swap_endian");
		
		int levelmin = cf.getValue<unsigned>("setup","levelmin");
		int levelmax = cf.getValue<unsigned>("setup","levelmax");
        block_buf_size_ = pow(pow(2,levelmin),2); //Npage=nrow^2; Number of particles in each page
		
		YHe_ = cf.getValueSafe<double>("cosmology","YHe",0.248);
        gamma_ = cf.getValueSafe<double>("cosmology","gamma",5.0/3.0);
		// Set header
        std::string thead;
        thead=cf.getValueSafe<std::string>("output","header","ICs generated using MUSIC");
        strcpy(header_.head,thead.c_str()); // text for the header; any easy way to add also the version?
        std::string ws = " "; // Filling with blanks. Any better way?
        for (int i=thead.size(); i<45;i++) 
        {
            header_.head[i]=ws[0]; 
        }
        header_.aexpN = astart_;
        header_.aexp0 = header_.aexpN;
		header_.amplt = 0.0; // Amplitude of density fluctuations
		header_.astep = cf.getValue<double>("output","astep"); // Seems that this must also be in the config file 
		header_.istep = 0; // step (=0 in IC)
		header_.partw = 0.0; // mass of highest res particle. SEE BELOW
        	header_.TINTG = 0; //=0 in IC
        	header_.EKIN = 0.0; //SUM 0.5 * m_i*(v_i**2) in code units
        	header_.EKIN1 = 0; //=0 in IC
        	header_.EKIN2 = 0; //=0 in IC
        	header_.AU0 = 0; //=0 in IC
        	header_.AEU0 = 0; //=0 in IC
		header_.NROWC = pow(2,levelmin); // Number of particles in 1 dim (number of particles per page = NROW**2) 
		header_.NGRIDC = pow(2,levelmin); // Number of cells in 1 dim
        	header_.nspecies = 0; // number of dm species
		for( int ilevel=levelmax; ilevel>=(int)levelmin; --ilevel )
		{
        		header_.nspecies+=1; 
		}
		//header_.partw  SEE BELOW
			
	        header_.Nseed = 0; // random number used ( 0 for MUSIC? or set the random number used in the lowest level?)
        	header_.Om0 = cf.getValue<double>("cosmology","Omega_m"); //Omega_m
	        header_.Oml0 = cf.getValue<double>("cosmology","Omega_L"); //Omega_L
        	header_.hubble = cf.getValue<double>("cosmology","H0"); //hubble constant h=H/100
	        header_.Wp5 = 0.0; // 0.0
		header_.Ocurv = 1.0 - header_.Oml0 - header_.Om0; //
	        header_.Omb0 = cf.getValue<double>("cosmology","Omega_b");; // this parameter only appears in header in hydro runs
		for (int i=0;i<10;i++)
		{
			header_.wpart[i] = 0.0; // extras[0-9] part. masses from high res to low res (normalized to low res particle)
			header_.lpart[i] = 0; // extras[10-19] # particles from high res to low res cumulative!!! 
		}
		for (int i=0;i<header_.nspecies;i++)
		{
			if (!do_baryons_)
				header_.wpart[i] = 1.0/pow(8.0,(header_.nspecies-i-1)); //from high res to lo res // 8 should be changed for internal variable?
			else
				header_.wpart[i] = ((header_.Om0-omegab_)/header_.Om0)/pow(8.0,(header_.nspecies-i-1)); 
		}
		for (int i=0;i<80;i++)
		{
	        	header_.extras[i] = 0.0; //extras[20-99] 
		}
        header_.extras[13] = cf.getValueSafe<double>("cosmology","Omega_b",0.0); //->0 in IC; fix it?
        header_.extras[14] = cf.getValue<double>("cosmology","sigma_8"); //->0 in IC; fix it?
        header_.extras[15] = cf.getValue<double>("cosmology","nspec"); //->0 in IC; fix it? Slope of the Power spectrum
        header_.extras[79] = cf.getValue<double>("setup","boxlength"); 
		header_.partw = header_.wpart[0]; // mass of highest res particle. 8 should be changed for internal variable?

		ptf_.astep=header_.astep;
		LOGINFO("ART : done header info.");
        
	}


    
    void write_dm_mass( const grid_hierarchy& gh )
	{
		//... store all the meta data about the grid hierarchy in header variables
		npcdm_ = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		for (int i=0;i<header_.nspecies;i++)
		{
			header_.lpart[i] = gh.count_leaf_cells(gh.levelmax()-i, gh.levelmax()); //cumulative!!
		}
		
        //... write data for dark matter......
		size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		std::vector<T_store> temp_dat;
		temp_dat.reserve(block_buf_size_);
		
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_mass );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		
		size_t blksize = sizeof(T_store)*nptot;
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		size_t nwritten = 0;
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
		{
            double pmass = omegam_/(1ul<<(3*ilevel)); // this needs to be adjusted to have the right units
            
            if( do_baryons_ )
                pmass *= (omegam_-omegab_)/omegam_;
			
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
						{
							if( temp_dat.size() <  block_buf_size_ )
								temp_dat.push_back( pmass );	
							else
							{
								ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*block_buf_size_ );	
								nwritten += block_buf_size_;
								temp_dat.clear();
								temp_dat.push_back( pmass );	
							}
						}
		}
		
		if( temp_dat.size() > 0 )
		{	
			ofs_temp.write( (char*)&temp_dat[0], sizeof(T_store)*temp_dat.size() );		
			nwritten+=temp_dat.size();
		}
		
		if( nwritten != nptot )
			throw std::runtime_error("Internal consistency error while writing temporary file for DM masses");
		
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for DM masses");
        
        
        ofs_temp.close();
		LOGINFO("ART : done write dm masses info.");
    }
    
    void write_dm_position( int coord, const grid_hierarchy& gh )
	{
        size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
		
	    //coordinates are in the range 1 - (NGRID+1)
    	// so scale factor is  scaleX = Box/NGRID -> to Mpc/h (Box in Mpc/h) 
		double xfac = header_.NGRIDC; 

		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_pos+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		size_t blksize = sizeof(T_store)*nptot;
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		size_t nwritten = 0;
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
						{
							double xx[3];
							gh.cell_pos(ilevel, i, j, k, xx);
							
							//xx[coord] = fmod( (xx[coord]+(*gh.get_grid(ilevel))(i,j,k)) + 1.0, 1.0 ) - 0.5;
							xx[coord] = (xx[coord]+(*gh.get_grid(ilevel))(i,j,k))*xfac+1.0; 
							
							if( temp_data.size() < block_buf_size_ )
								temp_data.push_back( xx[coord] );
							else
							{
								ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
								nwritten += block_buf_size_;
								temp_data.clear();
								temp_data.push_back( xx[coord] );
							}
						}
		
		if( temp_data.size() > 0 )
		{	
			ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
			nwritten += temp_data.size();
		}
		
		if( nwritten != nptot )
			throw std::runtime_error("Internal consistency error while writing temporary file for positions");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for positions");
		
		ofs_temp.close();
    }
    
    void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
        size_t nptot = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		std::vector<T_store> temp_data;
		temp_data.reserve( block_buf_size_ );
		
                //In ART velocities are P = a_expansion*V_pec/(x_0H_0) 
		// where x_0 = comoving cell_size=Box/Ngrid;H_0 = Hubble at z=0
		// so scale factor to physical km/s is convV= BoxV/AEXPN/NGRID 
		// (BoxV is Box*100; aexpn=current expansion factor)
                //internal units of MUSIC: To km/s just multiply by Lbox
		double vfac =  (header_.aexpN*header_.NGRIDC)/(100.0);  
        
		char temp_fname[256];
		sprintf( temp_fname, "___ic_temp_%05d.bin", 100*id_dm_vel+coord );
		std::ofstream ofs_temp( temp_fname, std::ios::binary|std::ios::trunc );
		
		size_t blksize = sizeof(T_store)*nptot;
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		size_t nwritten = 0;
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( ! gh.is_refined(ilevel,i,j,k) )
						{
							if( temp_data.size() < block_buf_size_ )
								temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
							else 
							{
								ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*block_buf_size_ );
								nwritten += block_buf_size_;
								temp_data.clear();
								temp_data.push_back( (*gh.get_grid(ilevel))(i,j,k) * vfac );
							}
							
						}
		
		if( temp_data.size() > 0 )
		{	
			ofs_temp.write( (char*)&temp_data[0], sizeof(T_store)*temp_data.size() );
			nwritten += temp_data.size();
		}
		
		if( nwritten != nptot )
			throw std::runtime_error("Internal consistency error while writing temporary file for DM velocities");
		
		//... dump to temporary file
		ofs_temp.write( (char *)&blksize, sizeof(size_t) );
		
		if( ofs_temp.bad() )
			throw std::runtime_error("I/O error while writing temporary file for DM velocities");
		
		ofs_temp.close();
    }
    
    void write_dm_density( const grid_hierarchy& gh )
	{
		//... we don't care about DM density for art
	}
    
    void write_dm_potential( const grid_hierarchy& gh )
	{ }
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ }
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{
        
    }
	
    void write_gas_position( int coord, const grid_hierarchy& gh )
	{
        //... we don't care about gas positions in art
    }
    
    void write_gas_density( const grid_hierarchy& gh )
	{
        
    }

	void finalize( void )
	{ 
        this->write_header_file();
        this->write_pt_file();    
		this->assemble_DM_file();
	}
};

namespace{
	output_plugin_creator_concrete<art_output_plugin<float> > creator("art");
}
