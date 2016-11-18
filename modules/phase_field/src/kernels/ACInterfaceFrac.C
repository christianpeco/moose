/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ACInterfaceFrac.h"

template<>
InputParameters validParams<ACInterfaceFrac>()

{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the kernel");
  params.addParam<MaterialPropertyName>("toughness_name", "Gc", "The fracture toughness used with the kernel");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  params.addParam<bool>("variable_L", true, "The mobility is a function of any MOOSE variable (if this is set to false L must be constant over the entire domain!)");
   params.addParam<bool>("variable_Gc", true, "The toughness is a function of any MOOSE variable (if this is set to false Gc must be constant over the entire domain!)");
  return params;
}

ACInterfaceFrac::ACInterfaceFrac(const InputParameters & parameters) :
    DerivativeMaterialInterface<JvarMapKernelInterface<Kernel> >(parameters),
    _L(getMaterialProperty<Real>("mob_name")),
    _kappa(getMaterialProperty<Real>("kappa_name")),
    _Gc(getMaterialProperty<Real>("toughness_name")),
    _variable_L(getParam<bool>("variable_L")),
    _variable_Gc(getParam<bool>("variable_Gc")),
    _dLdop(getMaterialPropertyDerivative<Real>("mob_name", _var.name())),
    _d2Ldop2(getMaterialPropertyDerivative<Real>("mob_name", _var.name(), _var.name())),
    _dkappadop(getMaterialPropertyDerivative<Real>("kappa_name", _var.name())),
    _dGcdop(getMaterialPropertyDerivative<Real>("toughness_name", _var.name())),
    _d2Gcdop2(getMaterialPropertyDerivative<Real>("toughness_name", _var.name(), _var.name())),
    _nvar(_coupled_moose_vars.size()),
    _dLdarg(_nvar),
    _d2Ldargdop(_nvar),
    _d2Ldarg2(_nvar),
    _dkappadarg(_nvar),
    _dGcdarg(_nvar),
    _d2Gcdargdop(_nvar),
    _d2Gcdarg2(_nvar),
    _gradarg(_nvar)
{
  // Get mobility and kappa derivatives and coupled variable gradients
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariable * ivar = _coupled_moose_vars[i];
    const VariableName iname = ivar->name();
    if (iname == _var.name())
      mooseError("The kernel variable should not be specified in the coupled `args` parameter.");

    _dLdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname);
    _dkappadarg[i] = &getMaterialPropertyDerivative<Real>("kappa_name", iname);
    _dGcdarg[i] = &getMaterialPropertyDerivative<Real>("toughness_name", iname);

    _d2Ldargdop[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _var.name());
    _d2Gcdargdop[i] = &getMaterialPropertyDerivative<Real>("toughness_name", iname, _var.name());

    _gradarg[i] = &(ivar->gradSln());

    _d2Ldarg2[i].resize(_nvar);
    for (unsigned int j = 0; j < _nvar; ++j)
      _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _coupled_moose_vars[j]->name());

    _d2Gcdarg2[i].resize(_nvar);
    for (unsigned int j = 0; j < _nvar; ++j)
      _d2Gcdarg2[i][j] = &getMaterialPropertyDerivative<Real>("toughness_name", iname, _coupled_moose_vars[j]->name());
  }
}

void
ACInterfaceFrac::initialSetup()
{
  validateCoupling<Real>("mob_name");
  validateCoupling<Real>("kappa_name");
  validateCoupling<Real>("toughness_name");
}

RealGradient
ACInterfaceFrac::gradL()
{
  RealGradient g = _grad_u[_qp] * _dLdop[_qp];
  for (unsigned int i = 0; i < _nvar; ++i)
    g += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];
  return g;
}

RealGradient
ACInterfaceFrac::gradGc()
{
  RealGradient g = _grad_u[_qp] * _dGcdop[_qp];
  for (unsigned int i = 0; i < _nvar; ++i)
    g += (*_gradarg[i])[_qp] * (*_dGcdarg[i])[_qp];
  return g;
}

RealGradient
ACInterfaceFrac::kappaNablaLGcPsi()
{
  // sum is the product rule gradient \f$ \nabla (L\psi) \f$
  RealGradient sum = _L[_qp] * _Gc[_qp] *_grad_test[_i][_qp];

  if (_variable_L)
    sum += gradL() * _Gc[_qp] * _test[_i][_qp];

  if (_variable_Gc)
    sum += _L[_qp] * gradGc() * _test[_i][_qp];

  return _kappa[_qp] * sum;
}

Real
ACInterfaceFrac::computeQpResidual()
{
  return _grad_u[_qp] * kappaNablaLGcPsi();
}

Real
ACInterfaceFrac::computeQpJacobian()
{
  // dsum is the derivative \f$ \frac\partial{\partial \eta} \left( \nabla (L\psi) \right) \f$
  RealGradient dsum = (_dkappadop[_qp] * _L[_qp] * _Gc[_qp] + _kappa[_qp] * _dLdop[_qp] * _Gc[_qp] + _kappa[_qp] * _L[_qp] * _dGcdop[_qp]) * _phi[_j][_qp] * _grad_test[_i][_qp];

  // compute the derivative of the gradient of the mobility
  if (_variable_L)
  {
    RealGradient dgradL =   _grad_phi[_j][_qp] * _dLdop[_qp]
                          + _grad_u[_qp] * _phi[_j][_qp] * _d2Ldop2[_qp];

    for (unsigned int i = 0; i < _nvar; ++i)
      dgradL += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Ldargdop[i])[_qp];

    dsum += (_kappa[_qp] * dgradL * _Gc[_qp] + _kappa[_qp] * gradL() * _dGcdop[_qp] + _dkappadop[_qp] * _phi[_j][_qp] * gradL() * _Gc[_qp]) * _test[_i][_qp];
  }

  if (_variable_Gc)
  {
    RealGradient dgradGc =   _grad_phi[_j][_qp] * _dGcdop[_qp]
                          + _grad_u[_qp] * _phi[_j][_qp] * _d2Gcdop2[_qp];

    for (unsigned int i = 0; i < _nvar; ++i)
      dgradGc += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Gcdargdop[i])[_qp];

    dsum += (_kappa[_qp] * _dLdop[_qp] * gradGc() + _kappa[_qp] * _L[_qp] * dgradGc + _dkappadop[_qp] * _phi[_j][_qp] * _L[_qp] * gradGc()) * _test[_i][_qp];
  }

  return _grad_phi[_j][_qp] * kappaNablaLGcPsi() + _grad_u[_qp] * dsum;
}

Real
ACInterfaceFrac::computeQpOffDiagJacobian(unsigned int jvar)
{
  // get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);

  // dsum is the derivative \f$ \frac\partial{\partial \eta} \left( \nabla (L\psi) \right) \f$
  RealGradient dsum = ((*_dkappadarg[cvar])[_qp] * _L[_qp] * _Gc[_qp] + _kappa[_qp] * (*_dLdarg[cvar])[_qp] * _Gc[_qp] + _kappa[_qp] * _L[_qp] * (*_dGcdarg[cvar])[_qp] ) * _phi[_j][_qp] * _grad_test[_i][_qp];

  // compute the derivative of the gradient of the mobility
  if (_variable_L)
  {
    RealGradient dgradL =   _grad_phi[_j][_qp] * (*_dLdarg[cvar])[_qp]
                          + _grad_u[_qp] * _phi[_j][_qp] * (*_d2Ldargdop[cvar])[_qp];

    for (unsigned int i = 0; i < _nvar; ++i)
      dgradL += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Ldarg2[cvar][i])[_qp];

    dsum += (_kappa[_qp] * dgradL * _Gc[_qp] + _kappa[_qp] * gradL() * _dGcdop[_qp] + _dkappadop[_qp] * _phi[_j][_qp] * gradL() * _Gc[_qp]) * _test[_i][_qp];
  }

  if (_variable_Gc)
  {
    RealGradient dgradGc =   _grad_phi[_j][_qp] * (*_dGcdarg[cvar])[_qp]
                          + _grad_u[_qp] * _phi[_j][_qp] * (*_d2Gcdargdop[cvar])[_qp];

    for (unsigned int i = 0; i < _nvar; ++i)
      dgradGc += (*_gradarg[i])[_qp] * _phi[_j][_qp] * (*_d2Gcdarg2[cvar][i])[_qp];

    dsum += (_kappa[_qp] * _dLdop[_qp] * gradGc()  + _kappa[_qp] * _L[_qp] * dgradGc + _dkappadop[_qp] * _phi[_j][_qp] *_L[_qp] * gradGc()) * _test[_i][_qp];
  }


  return _grad_u[_qp] * dsum;
}
