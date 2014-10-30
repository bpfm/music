/*
 * output_bigfile.cc - This file is part of MUSIC -
 * a code to generate multi-scale initial conditions 
 * for cosmological simulations
 * 
 * Copyright (C) 2010  Oliver Hahn
 * 
 * Plugin: Yu Feng(yfeng1@berkeley.edu)
 *
 * ( Modified from output_arepo.cc )
 */

#ifdef HAVE_BIGFILE
 
#define GAS_PARTTYPE 0
#define HIGHRES_DM_PARTTYPE 1
#define COARSE_DM_DEFAULT_PARTTYPE 2
#define STAR_PARTTYPE 4
#define NTYPES 6

#include <sstream>
#include <string>
#include <algorithm>
#include "output.hh"

extern "C" {
#include <bigfile.h>
}

template<typename T>
const char * typestring( void )
{
  if( typeid(T) == typeid(int) )
    return "i4";

  if( typeid(T) == typeid(unsigned) )
    return "u4";

  if( typeid(T) == typeid(float) )
    return "f4";

  if( typeid(T) == typeid(double) )
    return "f8";
  
  if( typeid(T) == typeid(long long) )
		return "i8";
	
  if( typeid(T) == typeid(unsigned long long) )
		return "u8";
	
  if( typeid(T) == typeid(size_t) )
		return "u8";
  
  std::cerr << " - Error: [HDF_IO] trying to evaluate unsupported type in typestring\n\n";
  return "";
}
class bigfile_output_plugin : public output_plugin
{ 
protected:
	
	// header/config
	std::vector<int> nPart;
	std::vector<long long> nPartTotal;
	std::vector<double> massTable;
	double time, redshift, boxSize;
	int numFiles, doublePrec;
	
	double omega0, omega_L, hubbleParam;
	
	// configuration
	double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
	double omega_b, rhoCrit;
	double posFac, velFac;
	int coarsePartType, nPartTotAllTypes;
	bool doBaryons, useLongIDs;
	
	size_t npfine, npart, npcoarse;
	std::vector<size_t> levelcounts;
	
	// parameter file hints
	int pmgrid, gridboost;
	double softening, Tini;
	
	using output_plugin::cf_;
	
    template <typename T>
	void bigfilehelper(std::string& fieldName, int partTypeNum, const std::vector<T> &data, int nmemb, bool write ) {
        char * bname = (char*) alloca(128);
        sprintf(bname, "%d/%s", partTypeNum, fieldName.c_str()); 
        BigFile bf = {0}; 
        BigBlock block = {0}; 
        BigBlockPtr ptr= {0}; 
        BigArray array = {0}; 

        size_t dims[] = {data.size() / nmemb, nmemb};
        size_t * fsize = (size_t *) alloca(sizeof(size_t) * numFiles);
        size_t total = data.size() / nmemb;
        int i;
        for (i = 0; i < numFiles; i ++) {
            fsize[i] = (i + 1) * total / numFiles - (i * total / numFiles);
        }

        big_array_init(&array, (void *)&data[0], typestring<T>(), 2, dims, NULL);
        big_block_seek(&block, &ptr, 0);
        big_file_open(&bf, fname_.c_str());
        if (write) {
            big_file_create_block(&bf, &block, bname, typestring<T>(), nmemb, numFiles, fsize);
            big_block_write(&block, &ptr, &array);
        } else {
            big_file_open_block(&bf, &block, bname);
            big_block_read(&block, &ptr, &array);
        }
        big_block_close(&block);
        big_file_close(&bf);
    }

	// Nx1 vector (e.g. masses,particleids)
	template< typename T >
	void writeHDF5_a( std::string fieldName, int partTypeNum, const std::vector<T> &data )
	{
        bigfilehelper(fieldName, partTypeNum, data, 1, true);
	}
	
	// Nx3 vector (e.g. positions,velocities), where coord = index of the second dimension (writen one at a time)
    template<typename T>
	void writeHDF5_b( std::string fieldName, int coord, int partTypeNum, std::vector<T> &data, bool readFlag = false )
	{
        bigfilehelper(fieldName, partTypeNum, data, 3, !readFlag);
	}
	
	// called from finalize()
  void generateAndWriteIDs( void )
  {
		long long offset = 0;
		nPartTotAllTypes = 0;
		
		for( size_t i=0; i < nPartTotal.size(); i++ )
		{
			if( !nPartTotal[i] )
				continue;
				
			nPartTotAllTypes += nPartTotal[i];
				
		    std::vector<long long> ids = std::vector<long long>(nPartTotal[i]);

            for( long long j=0; j < nPartTotal[i]; j++ )
                ids[j] = offset + j;
                
            writeHDF5_a( "ID", i, ids );
        
			// make IDs of all particle types sequential (unique) = unnecessary, but consistent with gadget output format
			offset += nPartTotal[i];
		}
	}
	
	void countLeafCells( const grid_hierarchy& gh )
	{
		npfine = 0; npart = 0; npcoarse = 0;
		
		npfine = gh.count_leaf_cells(gh.levelmax(), gh.levelmax());
		npart = gh.count_leaf_cells(gh.levelmin(), gh.levelmax());
		
		if( levelmax_ != levelmin_ ) // multimass
			npcoarse = gh.count_leaf_cells(gh.levelmin(), gh.levelmax()-1);
	}

public:
	bigfile_output_plugin( config_file& cf ) : output_plugin( cf )
	{
		// ensure that everyone knows we want to do SPH, implies: bsph=1, bbshift=1, decic_baryons=1
		// -> instead of just writing gas densities (which are here ignored), the gas displacements are also written
		cf.insertValue("setup","do_SPH","yes");
		
		// init header and config parameters
		nPart      = std::vector<int>(NTYPES,0);
		nPartTotal = std::vector<long long>(NTYPES,0);
		massTable  = std::vector<double>(NTYPES,0.0);
		
		coarsePartType   = cf.getValueSafe<unsigned>("output","bigfile_coarsetype",COARSE_DM_DEFAULT_PARTTYPE);
		UnitLength_in_cm = cf.getValueSafe<double>("output","bigfile_unitlength",3.085678e21); // 1.0 kpc
		UnitMass_in_g    = cf.getValueSafe<double>("output","bigfile_unitmass",1.989e43); // 1.0e10 solar masses
		UnitVelocity_in_cm_per_s = cf.getValueSafe<double>("output","bigfile_unitvel",1e5); // 1 km/sec
		
	  omega0     = cf.getValue<double>("cosmology","Omega_m");
		omega_b    = cf.getValue<double>("cosmology","Omega_b");
		omega_L    = cf.getValue<double>("cosmology","Omega_L");
		redshift   = cf.getValue<double>("setup","zstart");
		boxSize    = cf.getValue<double>("setup","boxlength");
		doBaryons  = cf.getValueSafe<bool>("setup","baryons",false);
		numFiles   = cf.getValueSafe<unsigned>("output","bigfile_num_files",1);
		
		
		// factors which multiply positions and velocities
		time   = 1.0/(1.0+redshift);
		posFac = 3.085678e24 / UnitLength_in_cm; // MUSIC uses Mpc internally, i.e. posFac=1e3 for kpc output
		velFac = ( 1.0f / sqrt(time) ) * boxSize; // TODO: should be normalized by posFac?
		
		// critical density
		rhoCrit = 27.7519737e-9; // in h^2 1e10 M_sol / kpc^3
		rhoCrit *= pow(UnitLength_in_cm/3.085678e21, 3.0);
		rhoCrit *= (1.989e43/UnitMass_in_g);
		
		// calculate PMGRID suggestion
		pmgrid = pow(2,levelmin_) * 2; // unigrid
		gridboost = 1;
		
		if( levelmin_ != levelmax_ )
		{
			double lxref[3], x0ref[3], x1ref[3];
			double pmgrid_new;
			
			the_region_generator->get_AABB(x0ref,x1ref,levelmax_); // generalized beyond box
			for (int i=0; i < 3; i++)
			  lxref[i] = x1ref[i] - x0ref[i];
			
			// fraction box length of the zoom region
			lxref[0] = pow( (lxref[0]*lxref[1]*lxref[2]),0.333 );
			
			pmgrid_new = pow(2,levelmax_) * 2; // to cover entire box at highest resolution
			pmgrid_new *= lxref[0]; // only need to cover a fraction
			
			if( (gridboost=round(pmgrid_new/pmgrid)) > 1 )
				gridboost = pow(2, ceil(log(gridboost)/log(2.0))); // round to nearest, higher power of 2
		}
		
		// calculate Tini for gas
		hubbleParam = cf.getValue<double>("cosmology","H0")/100.0;
		
		double astart = 1.0/(1.0+redshift);
		double h2     = hubbleParam*hubbleParam;
		double adec   = 1.0/( 160.0*pow(omega_b*h2/0.022,2.0/5.0) );
		double Tcmb0  = 2.726;
		
		Tini = astart<adec? Tcmb0/astart : Tcmb0/astart/astart*adec;
		
		// calculate softening suggestion
		softening = (boxSize * posFac) / pow(2,levelmax_) / 40.0;
		
		// header and sanity checks
		if ( !doBaryons )
			massTable[HIGHRES_DM_PARTTYPE] = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmax_);
		else
			massTable[HIGHRES_DM_PARTTYPE] = (omega0-omega_b) * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmax_);
		
		if ( coarsePartType == GAS_PARTTYPE || coarsePartType == HIGHRES_DM_PARTTYPE)
      throw std::runtime_error("Error: Specified illegal Arepo particle type for coarse particles.");
		if ( coarsePartType == STAR_PARTTYPE )
			LOGWARN("WARNING: Specified coarse particle type will collide with stars if USE_SFR enabled.");

        BigFile bf = {0};
        big_file_create(&bf, fname_.c_str());    
	    big_file_close(&bf);	
	}
	
	~bigfile_output_plugin()
	{	}
	
	/* ------------------------------------------------------------------------------- */
	
	void write_dm_mass( const grid_hierarchy& gh )
	{
		countLeafCells(gh);
		
		// fill levelcount for header
		levelcounts = std::vector<size_t>(levelmax_-levelmin_+1);
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			levelcounts[gh.levelmax()-ilevel] = gh.count_leaf_cells(ilevel, ilevel);
		
    if( levelmax_ > levelmin_ +1 ) // morethan2bnd
		{
			// DM particles will have variable masses
			size_t count = 0;
			
			std::vector<float> data(npcoarse);
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
			{
        // baryon particles live only on finest grid, these particles here are total matter particles
				float pmass = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*ilevel);	
				
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
							{
								data[count++] = pmass;
							}
			}
			
			if( count != npcoarse )
				throw std::runtime_error("Internal consistency error while writing masses");
				
			writeHDF5_a( "Mass", coarsePartType, data ); // write DM
			
		}
		else
		{
			// DM particles will all have the same mass, just write to massTable
			if( levelmax_ != levelmin_ ) // multimass
			  massTable[coarsePartType] = omega0 * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmin_);
		}		
	}
	
	void write_dm_position( int coord, const grid_hierarchy& gh )
	{
		countLeafCells(gh);
		
		// update header
		nPart[HIGHRES_DM_PARTTYPE] = npfine;
		nPart[coarsePartType]      = npcoarse;
		nPartTotal[HIGHRES_DM_PARTTYPE] = npfine;
		nPartTotal[coarsePartType]      = npcoarse;
		
		// FINE: collect displacements and convert to absolute coordinates with correct units
		int ilevel = gh.levelmax();
		
		std::vector<float> data(npfine);
		size_t count = 0;
		
		for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
			for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
				for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
					if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
					{
						double xx[3];
						gh.cell_pos(ilevel, i, j, k, xx);
							
						xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
						xx[coord] = fmod( xx[coord] + boxSize,boxSize );
						
						data[count++] = (float) (xx[coord] * posFac);
					}
						
		writeHDF5_b( "Position", coord, HIGHRES_DM_PARTTYPE, data );	// write fine DM
		
		if( count != npfine )
			throw std::runtime_error("Internal consistency error while writing fine DM pos");
		
		// COARSE: collect displacements and convert to absolute coordinates with correct units
		if( levelmax_ != levelmin_ ) // multimass
		{
			data = std::vector<float> (npcoarse,0.0);
			count = 0;
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
							{
								double xx[3];
								gh.cell_pos(ilevel, i, j, k, xx);
								
								xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
								
								if ( !doBaryons ) // if so, we will handle the mod in write_gas_position
									xx[coord] = fmod( xx[coord] + boxSize,boxSize ) * posFac;
																
								data[count++] = (float) xx[coord];
							}
				
				if( count != npcoarse )
					throw std::runtime_error("Internal consistency error while writing coarse DM pos");
					
				writeHDF5_b( "Position", coord, coarsePartType, data ); // write coarse DM
		}
	}
	
	void write_dm_velocity( int coord, const grid_hierarchy& gh )
	{
		countLeafCells(gh);
			
		// FINE: collect velocities and convert to correct units
		int ilevel = gh.levelmax();
		
		std::vector<float> data(npfine);
		size_t count = 0;
		
		for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
			for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
				for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
					if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
					{
						data[count++] = (*gh.get_grid(ilevel))(i,j,k) * velFac;
					}
						
		writeHDF5_b( "Velocity", coord, HIGHRES_DM_PARTTYPE, data ); // write fine DM
		
		if( count != npfine )
			throw std::runtime_error("Internal consistency error while writing fine DM pos");
		
		// COARSE: collect velocities and convert to correct units
		if( levelmax_ != levelmin_ ) // multimass
		{
			data = std::vector<float> (npcoarse,0.0);
			count = 0;
			
			for( int ilevel=gh.levelmax()-1; ilevel>=(int)gh.levelmin(); --ilevel )
				for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
					for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
						for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
							if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
							{
								data[count++] = (*gh.get_grid(ilevel))(i,j,k) * velFac;
							}
				
				if( count != npcoarse )
					throw std::runtime_error("Internal consistency error while writing coarse DM pos");
					
				writeHDF5_b( "Velocity", coord, coarsePartType, data ); // write coarse DM
		}
	
	}
	
	void write_dm_density( const grid_hierarchy& gh )
	{ /* skip */ }
	
	void write_dm_potential( const grid_hierarchy& gh )
	{ /* skip */ }
	
	/* ------------------------------------------------------------------------------- */
	
	void write_gas_velocity( int coord, const grid_hierarchy& gh )
	{	
		countLeafCells(gh);
		
		std::vector<float> gas_data(npart); // read/write gas at all levels from the gh
		size_t count = 0;
		
		for( int ilevel=levelmax_; ilevel>=(int)levelmin_; --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
						{
							gas_data[count++] = (*gh.get_grid(ilevel))(i,j,k) * velFac;
						}
						
		if( count != npart )
			throw std::runtime_error("Internal consistency error while writing GAS pos");
					
		// calculate modified DM velocities if: multimass and baryons present
		if( doBaryons && npcoarse )
		{
			double facb = omega_b / omega0;
			double facc = (omega0 - omega_b) / omega0;
			
			std::vector<float> dm_data(npcoarse);
			
			writeHDF5_b( "Velocity", coord, coarsePartType, dm_data, true ); // read coarse DM vels
			
			// overwrite 
			for( size_t i=0; i < npcoarse; i++ )
				dm_data[i] = facc*dm_data[i] + facb*gas_data[npfine + i];

			writeHDF5_b( "Velocity", coord, coarsePartType, dm_data ); // overwrite coarse DM vels
		} // dm_data deallocated
		
		// restrict gas_data to fine only and request write
		std::vector<float> data( gas_data.begin() + 0, gas_data.begin() + npfine );
		
		std::vector<float>().swap( gas_data ); // deallocate
		
		writeHDF5_b( "Velocity", coord, GAS_PARTTYPE, data );	 // write highres gas
	}
	
	void write_gas_position( int coord, const grid_hierarchy& gh )
	{
		countLeafCells(gh);
		
		// update header (will actually write only gas at levelmax)
		nPart[GAS_PARTTYPE] = npfine;
		nPartTotal[GAS_PARTTYPE] = npfine;
		
		std::vector<double> gas_data(npart); // read/write gas at all levels from the gh
		size_t count = 0;
		
		double h = 1.0/(1ul<<gh.levelmax());
		
		for( int ilevel=gh.levelmax(); ilevel>=(int)gh.levelmin(); --ilevel )
			for( unsigned i=0; i<gh.get_grid(ilevel)->size(0); ++i )
				for( unsigned j=0; j<gh.get_grid(ilevel)->size(1); ++j )
					for( unsigned k=0; k<gh.get_grid(ilevel)->size(2); ++k )
						if( gh.is_in_mask(ilevel,i,j,k) && !gh.is_refined(ilevel,i,j,k) )
						{
							double xx[3];
							gh.cell_pos(ilevel, i, j, k, xx);
							
							// shift particle positions (this has to be done as the same shift
							// is used when computing the convolution kernel for SPH baryons)
							xx[coord] += 0.5*h;
														
							xx[coord] = (xx[coord] + (*gh.get_grid(ilevel))(i,j,k)) * boxSize;
											
							gas_data[count++] = xx[coord];
						}
					
		if( count != npart )
			throw std::runtime_error("Internal consistency error while writing coarse DM pos");
					
		// calculate modified DM coordinates if: multimass and baryons present
		if( doBaryons && npcoarse )
		{
			double facb = omega_b / omega0;
			double facc = (omega0 - omega_b) / omega0;
			
			std::vector<float> dm_data(npcoarse);
			
			writeHDF5_b( "Position", coord, coarsePartType, dm_data, true ); // read coarse DM vels
			
			// overwrite 
			for( size_t i=0; i < npcoarse; i++ ) {
				dm_data[i] = facc*dm_data[i] + facb*gas_data[npfine + i];
				dm_data[i] = fmod( dm_data[i] + boxSize, boxSize ) * posFac;
			}

			writeHDF5_b( "Position", coord, coarsePartType, dm_data ); // overwrite coarse DM vels
		}
		
		// restrict gas_data to fine only and request write
		//std::vector<float> data( gas_data.begin() + 0, gas_data.begin() + npfine );
		
		std::vector<float> data(npfine);
		
		for( size_t i = 0; i < npfine; i++ )
			data[i] = (float) ( fmod( gas_data[i] + boxSize, boxSize ) * posFac );
		
		std::vector<double>().swap( gas_data ); // deallocate
		
		writeHDF5_b( "Position", coord, GAS_PARTTYPE, data ); // write highres gas

	}
	
	void write_gas_density( const grid_hierarchy& gh )
	{
		// if only saving highres gas, then all gas cells have the same initial mass
		// do not write out densities as we write out displacements
		if( doBaryons )
			massTable[GAS_PARTTYPE] = omega_b * rhoCrit * pow(boxSize*posFac,3.0)/pow(2,3*levelmax_);
	}
	
	void write_gas_potential( const grid_hierarchy& gh )
	{ /* skip */ }
	
	void finalize( void )
	{		
		// generate and add contiguous IDs for each particle type we have written
		generateAndWriteIDs();
		
		// write final header (some of these fields are required, others are extra info)
		BigFile bf = {0};
        BigBlock header = {0};
        big_file_open(&bf, fname_.c_str());    
        big_file_create_block(&bf, &header, "header", NULL, 0, 0, NULL);

		big_block_set_attr(&header, "TotNumPart",          &nPartTotal[0], "i8", 6);
		big_block_set_attr(&header, "MassTable",              &massTable[0], "f8", 6);
        double BoxSize = posFac * boxSize; /* convert to gadget unit*/
		big_block_set_attr(&header, "BoxSize",                &BoxSize, "f8", 1);
		big_block_set_attr(&header, "Time",                   &time, "f8", 1);
		big_block_set_attr(&header, "Redshift",               &redshift, "f8", 1);
		big_block_set_attr(&header, "Omega0",                 &omega0, "f8", 1);
		big_block_set_attr(&header, "OmegaLambda",            &omega_L, "f8", 1);
		big_block_set_attr(&header, "OmegaBaryon",            &omega_b, "f8", 1);
		big_block_set_attr(&header, "HubbleParam",            &hubbleParam, "f8", 1);
		big_block_set_attr(&header, "Music_levelmin",         &levelmin_, "i4", 1);
		big_block_set_attr(&header, "Music_levelmax",         &levelmax_, "i4", 1);
		big_block_set_attr(&header, "Music_levelcounts",      &levelcounts, "i8", 1);
		big_block_set_attr(&header, "suggested_pmgrid",       &pmgrid, "f8", 1);
		big_block_set_attr(&header, "suggested_gridboost",    &gridboost, "f8", 1);
		big_block_set_attr(&header, "suggested_highressoft",  &softening, "f8", 1);
		big_block_set_attr(&header, "suggested_gas_Tinit",    &Tini, "f8", 1);
		
        big_block_close(&header);
        big_file_close(&bf);
		// output particle counts
		std::cout << " - Arepo : wrote " << nPartTotAllTypes << " particles to file..." << std::endl;
		for( size_t i=0; i < nPartTotal.size(); i++ )
			std::cout << "    type [" << i << "] : " << std::setw(12) << nPartTotal[i] << std::endl;
			
		// give config/parameter file hints		
		if( NTYPES > 6 )
			std::cout << " - Arepo: Using [" << NTYPES << "] particle types, set NTYPES to match." << std::endl;
		if( doBaryons )
			std::cout << " - Arepo: Wrote gas, set REFINEMENT_HIGH_RES_GAS and GENERATE_GAS_IN_ICS with "
			          << "SPLIT_PARTICLE_TYPE=" << pow(2,coarsePartType) << "." << std::endl;
		if( levelmax_ != levelmin_ )
			std::cout << " - Arepo: Have zoom type ICs, set PLACEHIGHRESREGION=" << pow(2,HIGHRES_DM_PARTTYPE)
                << " (suggest PMGRID=" << pmgrid << " with GRIDBOOST=" << gridboost << ")." << std::endl;
		if( levelmax_ > levelmin_ + 1 )
			std::cout << " - Arepo: More than one coarse DM mass using same type, set INDIVIDUAL_GRAVITY_SOFTENING=" 
			          << pow(2,coarsePartType) << " (+" << pow(2,STAR_PARTTYPE) << " if including stars)." << std::endl;
		if( doBaryons )
			std::cout << " - Arepo: Set initial gas temperature to " << std::fixed << std::setprecision(3) << Tini << " K." << std::endl;
		std::cout << " - Arepo: Suggest grav softening = " << std::setprecision(3) << softening << " for high res DM." << std::endl;
			
	}
	
};

namespace{
	output_plugin_creator_concrete< bigfile_output_plugin > creator("bigfile");
}
#endif /*BIGFILE*/
