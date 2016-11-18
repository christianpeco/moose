/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SurfactantPressure2DFrac.h"
#include "Function.h"

template<>
InputParameters validParams<SurfactantPressure2DFrac>()
{
  InputParameters params = validParams<Material>();

  params.addParam<Real>("const", "Constant");
  params.addParam<int>("pres_type", "1-quadratic 0- linear");
  params.addParam<bool>("irrev", "Irreversibility of overstress 1-yes 0- no");
  params.addParam<bool>("staggered", "0-no 1-yes");

  params.addParam<FunctionName>("loadingfunction", "The loading function");
  params.addRequiredParam<std::string>("base_name", "Elastic source name for strain");
  params.addCoupledVar("concentration", "Concentration of substance generating the pressure");
  params.addCoupledVar("concentration_ini", "Initial profile generating the pressure");


  return params;
}

SurfactantPressure2DFrac::SurfactantPressure2DFrac(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),

    //Get parameters : viscosity, diameter, L_const

    _const_set( parameters.isParamValid("const") ),
    _pres_type_set( parameters.isParamValid("pres_type") ),
    _irrev_set( parameters.isParamValid("irrev") ),
    _staggered_set( parameters.isParamValid("staggered") ),

    _const( _const_set ? getParam<Real>("const") : 1.0 ),
    _pres_type( _pres_type_set ? getParam<int>("pres_type") : 0 ),
    _irrev( _irrev_set ? getParam<bool>("irrev") : true ),
    _staggered( _staggered_set ? getParam<bool>("staggered") : true ),

    _loadingfunction(getFunction("loadingfunction")),
    _base_name(getParam<std::string>("base_name")),


    // Get concentration and trace
    _concentration(coupledValue("concentration")),
    _concentration_ini(coupledValue("concentration_ini")),
    _trace(getMaterialPropertyByName<Real>(_base_name + "_" +  "trace_st")),
    _trace_old(getMaterialPropertyOldByName<Real>(_base_name + "_" + "trace_st")),

    _concentration_history(declareProperty<Real>(_base_name + "_" + "concentration_history")),
    _concentration_history_old(declarePropertyOld<Real>(_base_name + "_" + "concentration_history")),
    _trace_history(declareProperty<Real>(_base_name + "_" + "trace_history")),
    _trace_history_old(declarePropertyOld<Real>(_base_name + "_" + "trace_history")),


    //Declare overstress that kernels can use, and derivatives
    _overisostress(declareProperty<Real>(_base_name + "_" + "overisostress")),
    _energy_over(declareProperty<Real>(_base_name + "_" + "energy_over")),
    _d_energy_over_dc(declareProperty<Real>(_base_name + "_" + "d_energy_over_dc")),
    _d2_energy_over_dc2(declareProperty<Real>(_base_name + "_" + "d2_energy_over_dc2"))


{}


void
SurfactantPressure2DFrac::initQpStatefulProperties()
{
_trace_history[_qp]=0.0;
_concentration_history[_qp]=_concentration_ini[_qp];
}



void
SurfactantPressure2DFrac::computeQpProperties()
{

    if (_trace_history[_qp]<_trace_old[_qp]) _trace_history[_qp] = _trace_old[_qp];
    else _trace_history[_qp] = _trace_history_old[_qp];


 // Compute Marangoni flux forces (based on previous step always)
 // expansion stress (which can be irreversible o not)



    //create concentration history max
    if (_irrev==true)
    {

        Real topvalue = std::max(_concentration[_qp],_concentration_ini[_qp]);

        if (_concentration_history_old[_qp] <= topvalue)
        {
            _concentration_history[_qp] = topvalue;
        }

        else {_concentration_history[_qp] = _concentration_history_old[_qp];}
    }
    if (_irrev==false) {_concentration_history[_qp] = _concentration[_qp];}


    //cuts off

    if (_concentration_history[_qp]>1.0) _concentration_history[_qp]=1.0;
    if (_concentration_history[_qp]<0.0) _concentration_history[_qp]=0.0;



   //compute overstress (now staggered)
   if (_pres_type==0)//linear
   {
   _overisostress[_qp] =  - _loadingfunction.value(_t, _q_point[_qp])*_const*_concentration_history_old[_qp];

   _energy_over[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*(_concentration_history_old[_qp])*_trace_old[_qp];
   _d_energy_over_dc[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*_trace_old[_qp];
   _d2_energy_over_dc2[_qp] = 0.0;

   }
   else if (_pres_type==1)//quadratic
   {
   _overisostress[_qp] =  - _loadingfunction.value(_t, _q_point[_qp])*_const*_concentration_history_old[_qp]*_concentration_history_old[_qp]*0.5;

   _energy_over[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*(_concentration_history_old[_qp])*(_concentration_history_old[_qp])*0.5*_trace_old[_qp];
   _d_energy_over_dc[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*(_concentration_history_old[_qp])*_trace_old[_qp];
   _d2_energy_over_dc2[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*_trace_old[_qp];

   }
   else if (_pres_type==2)//interface (applied on the media, positive compressive stress)
   {
   _overisostress[_qp] =  + _loadingfunction.value(_t, _q_point[_qp])*_const*_concentration_history_old[_qp]*(1.0-_concentration_history_old[_qp]);

   _energy_over[_qp] = + _loadingfunction.value(_t, _q_point[_qp])*_const*_concentration_history_old[_qp]*(1.0 -_concentration_history_old[_qp])*_trace_old[_qp];
   _d_energy_over_dc[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*(1.0 -_concentration_history_old[_qp])*_trace_old[_qp];
   _d2_energy_over_dc2[_qp] = - _loadingfunction.value(_t, _q_point[_qp])*_const*_trace_old[_qp];

   }



}
