#ifndef MAXWELL_PARSER_H
#define MAXWELL_PARSER_H

#include <paramparser.h>
#include <maxwellparameters.h>

class MaxwellParser : ParamParser{
  public:
    MaxwellParser() : ParamParser(1){}

    virtual ~MaxwellParser(){}

    virtual Parameters& getParameters(std::string fname){
      MaxwellParameters params = MaxwellParameters();
      updateParameters(fname, &params);
      return params;
    }

    virtual void updateParameters(std::string fname, Parameters* params);
};

#endif
