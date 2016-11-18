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

#ifndef ReadIC_H
#define ReadIC_H

// MOOSE Includes
#include "InitialCondition.h"
#include <vector>

// Forward Declarations
class ReadIC;
template<>
InputParameters validParams<ReadIC>();


class ReadIC : public InitialCondition
{
public:
    ReadIC(const InputParameters & parameters);
    
    virtual Real value(const Point & p) override;
    virtual void initialSetup() override;
private:
    std::string _dataFile;
    std::vector<double> _vals;
    
};

#endif //ReadIC_H
