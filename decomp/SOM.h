#ifndef SOM_h__
#define SOM_h__
#include "SOMcomponent.h"
namespace DECOMP {

	const std::string VERSION = std::string("DECOMP++ compiled ") + std::string(__DATE__) + " - " + std::string(__TIME__);

	/// @brief A class representing Soil Organic Matter (SOM) with the decomposition properties from Wallman 2006 (https://doi.org/10.1016/j.envsoft.2004.09.026)
	///
	/// SOM consists of different C-pools (SOMcomponent) and a single N-pool.
	/// Decomposition of SOM is calculated for each pool. The net decomposition rate is used to
	/// calculate the N mineralisation rate.
	/// 
	/// The carbon balance of each SOM component i is given as:
	/// \f[ \frac{dC_i}{dt} = C_{in,i} - r_i(T, \theta, pH) C_i \f]
	///
	/// Where \f$\frac{dC_i}{dt}\f$ is the change rate of the Carbon stock of component i \f$C_i\f$ in kg/day. 
	///        \f$C_{in,i}\f$ is the added or created Carbon in flux in kg/day and \f$r_i\f$ is the decay rate of the component,
	///		   a function of soil temperature \f$T\f$, water content \f$\theta\f$ and soil water \f$pH\f$, given in Wallman 2006, eq. 10-13
	///
	/// From the decomposition rate of each component, we can calculate the total mineralisation rate as the mass weigthed mean 
	/// of the decomposition rate of each component:
	///
	/// \f[r_{tot}=\frac {\sum_{i=0}^n{r_i C_i}} {\sum_{i=0}^{n}{C_i}} \f]
	///
	/// The **Nitrogen** balance \f$N\f$ for the whole SOM is calculated from the effective mineralisation rate \f$r_{tot}\f$ over
	/// all components and a function of the C/N ratio of the soil, that indicates an N immobilising soil with a wide C/N ratio
	/// or a N releasing soil with a more narrow C/N ratio.
	///
	/// \f[\frac{dN}{dt} = N_{in} - r_{tot} N + r_{tot} N \frac{CN - CN_{min}}{CN_{max} - CN_{min}} \f]
	/// 
	/// Where \f$N_{in}\f$ is the input rate, CN the actual C/N ratio of the SOM and \f$CN_{min}, CN_{max}\f$ is the range of
	/// C/N ratios between the system shifts from N aquiring to an N releasing system.
	///
	class SOM
	{
	private:
		std::valarray<double> C_pools;
		static component_set pool_types;
	public:
		static component_set& get_pool_types();
		static SOMcomponent add_component(std::string name, bool is_stored,double k_pot, double E_a, double K_w, double n_w, double K_pH, double m_pH);

		double 
			CNmin, ///< Minimal natural C/N ratio (default 15) (needed for N immobilisation)
			CNmax; ///< Maximum natural C/N ration (default 40) (needed for N immobilisation)
		/// N-content
		double N;
		double get_C_pool(int index) const;
		void set_C_pool(int index, double pool_size);

#ifndef SWIG
		double& operator[](const SOMcomponent& component)
		{
			return this->C_pools[component.Id];
		}

		SOM& operator=(const SOM& copy);
#endif
		/// Returns the stored carbon of a single component
		/// Returns the sum of the C-pools
		double get_C_pool() const;
		/// Returns the C/N ratio
		double get_CN() const;

		SOM& operator*=(double right);
		SOM operator*(double right) const;		 
		SOM& operator+=(const SOM& right);
		SOM& operator-=(const SOM& right);
		SOM operator+(const SOM& right) const
		{
			SOM res= *this;
			res += right;
			return res;
		}
		SOM operator-(const SOM& right) const
		{
			SOM res= *this;
			res -= right;
			return res;
		}
		SOM& operator/=(double right);
		SOM operator/(double right)
		{
			SOM res= *this;
			res /= right;
			return res;
		}
		

		/// Returns the change rate of the pools
		/// @returns The change rate as an SOM object
		/// @param T Temperature in °C
		/// @param wetness Wetness in m3/m3
		/// @param pH pH-Value of the soil
		/// @param Nsol reactive N in the soil solution in kg
		SOM dCdt(double T, double wetness, double pH,double Nsol=0) const;

		
		SOM(const SOM& copy);
				
		/// Creates a new organic matter object
		/// @param N Nitrogen content (in mass)
		/// @param EDC Easily decomposable components (in mass)
		/// @param CELL Cellulose and similar components (in mass)
		/// @param LIGN Lignin and similar components (in mass)
		/// @param RC Resistant components (in mass)
		/// @param DOC dissolved components (in mass)
		SOM(double N=0.0, double EDC=0.0,double CELL=0.0, double LIGN=0.0, double RC=0.0, double DOC=0.0);

		
		SOM integrate(double dt, double T, double wetness, double pH);

		std::string to_string() const;
	};
	SOM operator*(double left, const SOM& right);

	SOM wood_litter();
	SOM leave_litter();
	SOM root_litter();
	SOM pure_DOC();



}

#endif // SOM_h__
