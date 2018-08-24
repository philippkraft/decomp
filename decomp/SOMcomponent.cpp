/*
 * SOMcomponent.cpp
 *
 *  Created on: 21.09.2009
 *      Author: philkraf
 */

#include "SOMcomponent.h"
#include <cmath>
namespace DECOMP {

  double SOMcomponent::decomp( double T, double wetness, double pH)	const
    {
        return k_pot/365.25 * f_Temp(T) * f_wet(wetness) * f_pH(pH);
    }

  SOMcomponent::SOMcomponent( std::string name,bool _is_stored,
			double _k_pot,double _E_a,
			double _K_w, double _n_w,
			double _K_pH, double _m_pH)
	: Id(count++), Name(name), k_pot(_k_pot), E_a(_E_a),
	K_w(_K_w), n_w(_n_w), K_pH(_K_pH), m_pH(_m_pH), is_stored(_is_stored)
	{
	}

  double SOMcomponent::f_Temp(double T) const
	{
		double
			R=8.314*0.001,
			T_R = 5.0,
			arr_gamma=this->E_a/(R*(T_R+273.16))-this->E_a/(R*(T+273.16));
		return exp(arr_gamma);
	}

  double SOMcomponent::f_pH(double pH) const
  {
  	double H_conc = pow(10.,-pH);
  	return 1.0/(1.0+K_pH * pow(H_conc,m_pH));

  }

  double SOMcomponent::f_wet(double wet) const
  {
  	return K_w*pow(wet,n_w)/(1.0+K_w*pow(wet,n_w));
  }

	SOMcomponent::SOMcomponent(const SOMcomponent & copy)
	: 	Id(copy.Id), Name(copy.Name),
		k_pot(copy.k_pot), E_a(copy.E_a),
		K_w(copy.K_w), n_w(copy.n_w), K_pH(copy.K_pH), m_pH(copy.m_pH),
		products(copy.products), is_stored(copy.is_stored)
    {

    }

  SOMcomponent::~SOMcomponent() {
  }
  component_set SOMcomponent::get_products() const
  {
  	component_set res;
  	for (product_map::const_iterator it=products.begin();it!=products.end();++it)
  	{
  		res.push_back(it->first);
  	}
  	return res;
  }

	int SOMcomponent::count(0);
}







