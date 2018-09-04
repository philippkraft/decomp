/*
 * SOMcomponent.h
 *
 *  Created on: 21.09.2009
 *      Author: philkraf
 */


#ifndef SOMCOMPONENT_H_
#define SOMCOMPONENT_H_

#include <string>
#include <map>
#include <vector>
#include <list>
#include <memory>
#include <valarray>

/// Set of components	


class SOM;
class SOMcomponent;
typedef std::vector<SOMcomponent> component_set;
/// Describes the properties of a component of the Soil Organic Matter
/// the original DECOMP model has 4 pool types: EDC, CELL, LIGN and RC and 2 flux types: CO2 and DOC
class SOMcomponent {
private:
    typedef std::map<SOMcomponent, double> product_map;
    product_map products;
    double f_Temp(double T) const;
    double f_wet(double wet) const;
    double f_pH(double pH) const;
public:
    SOMcomponent(std::string name, bool is_stored,
                double k_pot,double E_a,
                double K_w, double n_w,
                double K_pH, double m_pH=1.0);
    static int count;
    /// Id of the component
    const int Id;
    /// Name of the component
    std::string Name;
    /// True if the component is storable
    bool is_stored;

    double
        k_pot, ///< Potential decomposition rate of the component in 1/year
        E_a,   ///< Activation energy for the decomposition in kJ/mol
        K_w,   ///< Water function coefficient
        n_w,   ///< Water function exponent
        K_pH,  ///< Response coefficient in pH function in kmol/m3 (= mol/l)
        m_pH;  ///< Response exponent in pH function


    /// Set the fraction of a product
    void set_product(const SOMcomponent& product,double fraction);
    /// Get the fraction of a product
    double get_product_fraction(const SOMcomponent& product) const {
        product_map::const_iterator fract=products.find(product);
        if (fract!=this->products.end())
            return fract->second;
        else
            return 0.0;
    }
    /// Gets the list of products of this compound
    component_set get_products() const;
    /// Calculates the decomposition rate in 1/day of this component
    double decomp(double T, double wetness, double pH) const;
    SOMcomponent();
    bool operator<(const SOMcomponent& cmp) const { return Id<cmp.Id;}
    bool operator==(const SOMcomponent& cmp) { return Id==cmp.Id;}
    bool operator!=(const SOMcomponent& cmp) { return Id!=cmp.Id;}

#ifndef SWIG
    SOMcomponent& operator=(const SOMcomponent& copy);
 #endif
    /// Copy constructor
    SOMcomponent(const SOMcomponent& copy);
    ~SOMcomponent();
};



#endif /* SOMCOMPONENT_H_ */
