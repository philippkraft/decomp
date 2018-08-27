#include "SOM.h"
#include <algorithm>
#include <sstream>
#define min(a,b) ((a)<(b) ? (a) : (b))

SOM::SOM( const SOM& copy ) : C_pools(copy.C_pools),N(copy.N)	, CNmin(15.0), CNmax(40.0)
{
}

SOM& SOM::operator=( const SOM& copy )
{
    C_pools=copy.C_pools;
    N=copy.N;
    CNmin=copy.CNmin;
    CNmax=copy.CNmax;
    return *this;
}

SOM::SOM(double N, double EDC,double CELL, double LIGN, double RC, double DOC )
: N(N), CNmin(15.0), CNmax(40.0)
{
    C_pools.resize(SOMcomponent::count);
    C_pools[0]=EDC;
    C_pools[1]=LIGN;
    C_pools[2]=CELL;
    C_pools[3]=RC;
    C_pools[4]=DOC;
}

SOM SOM::dCdt( double T, double wetness, double pH, double Nsol ) const
{
    SOM dispatch;
    SOM decomp;
    const SOM& self = *this;
    for(component_set::const_iterator it = pool_types.begin(); it != pool_types.end(); ++it)
    {
        const SOMcomponent& comp=*it;

        double decomp_comp = self.get_C_pool(comp.Id) > 0 ? self.get_C_pool(comp.Id) * comp.decomp(T,wetness,pH) : 0.0;
        component_set products=comp.get_products();
        for(component_set::iterator p_it=products.begin();p_it!=products.end();++p_it)
        {
            dispatch[*p_it] += decomp_comp * comp.get_product_fraction(*p_it);
        }
        decomp[comp] = decomp_comp;
    }
    SOM result = dispatch - decomp;

    double C_pool=self.get_C_pool();

    if (C_pool>0 && N>0)
    {
        double
            net_min = -result.get_C_pool(),
            CN = self.get_CN(),
            grossNmin = net_min/CN,
            f_Nimmob=min(1,(CN - self.CNmin)/(self.CNmax - self.CNmin)),
            Nimmob = grossNmin * f_Nimmob;
        result.N = Nimmob - grossNmin;
    }
    else
        result.N =0.0;

    return result;
}


double SOM::get_C_pool() const
{
    double res=0.0;
    if (C_pools.size() != pool_types.size())
        throw std::runtime_error("DECOMP: Pool size array and pool type array out of sync!");
    for(component_set::const_iterator it = pool_types.begin(); it != pool_types.end(); ++it)
        if (it->is_stored) res+= C_pools[it->Id];
    return res;
}



double SOM::get_C_pool( int index ) const
{
    if (index>=0 && index < int(C_pools.size()))
        return C_pools[index];
    else
        throw std::out_of_range("DECOMP: Invalid component ID");

}
void SOM::set_C_pool( int index, double pool_size )
{
    if (index>=0 && index < int(C_pools.size()))
        C_pools[index] = pool_size;
    else
        throw std::out_of_range("DECOMP: Invalid component ID");

}

double SOM::get_CN() const
{
    return get_C_pool()/N;
}

SOM& SOM::operator*=( double right )
{
    C_pools *= right;
    N*=right;
    return *this;
}

SOM SOM::operator*( double right ) const
{
    SOM res=*this;
    res*=right;
    return res;
}

SOM& SOM::operator+=( const SOM& right )
{
    C_pools += right.C_pools;
    N  += right.N;
    return *this;
}

SOM& SOM::operator-=( const SOM& right )
{
    C_pools -= right.C_pools;
    N  -= right.N;
    return *this;
}

SOM& SOM::operator/=( double right )
{
    C_pools /= right;
    N  /= right;
    return *this;

}

SOM SOM::integrate( double dt, double T, double wetness, double pH )
{
    // Just a shortcut to this
    SOM& self=*this;

    // Calculate the change rate
    SOM rate = dCdt(T,wetness,pH);

    // Add the change rate to the current storage
    self += rate * dt;

    // Remove all components which are not storages from SOM
    for(component_set::const_iterator it = pool_types.begin(); it != pool_types.end(); ++it)
    {
        if (!it->is_stored)
            self[*it]=0.0;
    }

    // Prepare the calculated change rate for output components, by deleting stored compenents from rate
    for(component_set::const_iterator it = pool_types.begin(); it != pool_types.end(); ++it)
        if (it->is_stored)
            rate[*it]=0.0;
    // Adjust the sign for Nitrogen
    rate.N *= -1;
    return rate;
}

std::string SOM::to_string() const
{
    std::stringstream sstr;
    sstr.precision(4);
    sstr << "SOM(N=" << N;
    for(component_set::const_iterator it = SOM::pool_types.begin(); it != SOM::pool_types.end(); ++it)
    {
      double pool_size=get_C_pool(it->Id);
        if (pool_size>0)
            sstr << "," << it->Name << "=" << pool_size;
    }
    sstr << ")";
    return sstr.str();
}

const component_set& SOM::get_pool_types()
{
    return SOM::pool_types;
}

SOMcomponent SOM::add_component( std::string name, bool is_stored,double k_pot, double E_a, double K_w, double n_w, double K_pH, double m_pH )
{
    SOMcomponent new_comp(name,is_stored,k_pot,E_a,K_w,n_w,K_pH,m_pH);
    pool_types.push_back(new_comp);
    return new_comp;
}

component_set init_SOMcomponents()
{
    SOMcomponent
        EDC("EDC",true,240,18,9.4,3.4,65600),
        CELL("CELL",true,11,33,9.4,3.4,20500.),
        LIGN("LIGN",true,1.7,50,9.4,3.4,1050.),
        RC("RECALC",true,0.025,53,9.4,3.4,1050),
        DOC("DOC",false,0.025,50,110,2.454,20500),
        CO2("CO_2",false,0,0,0,0,0);


    EDC.set_product(CO2,0.45);
    CELL.set_product(CO2,0.45);
    LIGN.set_product(CO2,0.4);
    RC.set_product(CO2,0.4);
    DOC.set_product(CO2,0.75);

    EDC.set_product(DOC,0.45);
    CELL.set_product(DOC,0.45);
    LIGN.set_product(DOC,0.4);
    RC.set_product(DOC,0.5);

    LIGN.set_product(LIGN,0.1);

    EDC.set_product(RC,0.1);
    CELL.set_product(RC,0.1);
    LIGN.set_product(RC,0.1);
    RC.set_product(RC,0.1);
    DOC.set_product(RC,0.25);

    component_set res;
    res.push_back(EDC);
    res.push_back(CELL);
    res.push_back(LIGN);
    res.push_back(RC);
    res.push_back(DOC);
    res.push_back(CO2);
    return res;

}
component_set SOM::pool_types=init_SOMcomponents();

SOM wood_litter()
{
    return SOM(1./50., .04, .6, .27,.09);
}

SOM leave_litter()
{
    return SOM(1./50., .1,.5,.32,.08);
}

SOM root_litter( )
{
    return SOM(1./20., .21,.4,.33,.05);
}

SOM pure_DOC()
{
    return SOM(0.,0.,0.,0.,0.,1.0);
}

SOM operator*(double left, const SOM& right)
{
    SOM res=right;
    res *= left;
    return res;
}


