/*
 
 convolution_kernel.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
*/

#ifndef __CONVOLUTION_KERNELS_HH
#define __CONVOLUTION_KERNELS_HH

#include <string>
#include <map>

#include "config_file.hh"
#include "densities.hh"
#include "transfer_function.hh"


#define ACC_RF(i,j,k) (((((i)+nx)%nx)*ny+(((j)+ny)%ny))*2*(nz/2+1)+(((k)+nz)%nz))
#define ACC_RC(i,j,k) (((((i)+nxc)%nxc)*nyc+(((j)+nyc)%nyc))*2*(nzc/2+1)+(((k)+nzc)%nzc))

namespace convolution{

	//! encapsulates all parameters required for transfer function convolution
	struct parameters
	{
		int nx,ny,nz;
		double lx,ly,lz;//,boxlength;
		config_file *pcf;
		transfer_function* ptf;
		unsigned coarse_fact;
		bool deconvolve;
		bool is_finest;
		bool smooth;
	};
	
	
	/////////////////////////////////////////////////////////////////
	
	
	//! abstract base class for a transfer function convolution kernel
	class kernel{
	public:
		
		//! all parameters (physical/numerical)
		parameters cparam_;
		
		config_file *pcf_;
		transfer_function* ptf_;
		refinement_hierarchy* prefh_;
		tf_type type_;
		
		//! constructor
		kernel( config_file& cf, transfer_function* ptf, refinement_hierarchy& refh, tf_type type )
		: pcf_(&cf), ptf_(ptf), prefh_(&refh), type_(type)//cparam_( cp )
		{	}
		
		//! dummy constructor
		/*kernel( void )
		{	}*/
		
		//! compute/load the kernel
		virtual kernel* fetch_kernel( int ilevel, bool isolated=false ) = 0;
		
		//! virtual destructor
		virtual ~kernel(){ };
		
		//! purely virtual method to obtain a pointer to the underlying data
		virtual void* get_ptr() = 0;	
		
		//! free memory
		virtual void deallocate() = 0;
	};

	
	//! abstract factory class to create convolution kernels
	struct kernel_creator
	{
		//! creates a convolution kernel object
		virtual kernel * create( config_file& cf, transfer_function* ptf, refinement_hierarchy& refh, tf_type type ) const = 0;
		
		//! destructor
		virtual ~kernel_creator() { }
	};
	
	
	//! access map to the various kernel classes through the factory
	std::map< std::string, kernel_creator *>& get_kernel_map();
	
	
	//! actual implementation of the factory class for kernel objects
	template< class Derived >
	struct kernel_creator_concrete : public kernel_creator
	{
		//! constructor inserts the kernel class in the map
		kernel_creator_concrete( const std::string& kernel_name )
		{	get_kernel_map()[ kernel_name ] = this; 	}
		
		//! creates an instance of the kernel object
		kernel * create( config_file& cf, transfer_function* ptf, refinement_hierarchy& refh, tf_type type ) const
		{	return new Derived( cf, ptf, refh, type );	}
	};

	
	//! actual implementation of the FFT convolution (independent of the actual kernel)
	template< typename real_t >
	void perform( kernel* pk, void *pd );
	
	
	
	
} //namespace convolution

	
#endif //__CONVOLUTION_KERNELS_HH
