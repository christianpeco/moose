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


#include "ReadIC.h"

template<>
InputParameters validParams<ReadIC>()
{
    InputParameters params = validParams<InitialCondition>();
    params.addParam<std::string>("dataFile", "csv data file to be read in");
    params.addParam<std::vector<double>>("vals","Values read in from dataFile");
    return params;
}

ReadIC::ReadIC(const InputParameters & parameters) :
InitialCondition(parameters),
_dataFile(getParam<std::string>("dataFile")),
_vals(getParam<std::vector<double>>("vals"))
{}

void
ReadIC::initialSetup(){
    //read in data to get number of lines
    std::ifstream fileIN;
    
    fileIN.open(_dataFile);
    std::string lineA;
    int lines=0;
    while(fileIN.good()){
        while(getline(fileIN,lineA)){
            std::istringstream streamA(lineA);
            lines++;
        }
    }
    
    //string filename;
    fileIN.close();
    std::ifstream fileIN1;
    fileIN1.open(_dataFile);
    
    //
    //read data file
    std::vector<std::vector<double>> Values(lines,std::vector<double>(3));
    while(fileIN1.good()){
        
        for(int i=0;i<lines;i++){
            for (int j=0;j<3;j++){
                fileIN1 >> Values[i][j];
            }
        }
        
    }
    fileIN1.close();
    for(int i = 0; i < lines; i++){
        _vals.push_back(Values[i][2]);
        //std::cout<<_vals[i]<<std::endl;
    }
    

}
Real
ReadIC::value(const Point & p)
{

    return _vals[_current_node->id()];
    
}
