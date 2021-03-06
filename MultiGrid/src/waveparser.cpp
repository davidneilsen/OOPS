#include <waveparser.h>
#include <iostream>

using namespace std;

void WaveParser::updateParameters(string fname, Parameters* params){
  if(!checkId(params)){
    return;
  }
  reader.clearData();
  ParamReader::ParamResult result = reader.readFile(fname);
  if(result != ParamReader::SUCCESS){
    cout << "An error occurred while trying to read " << fname << ".\n";
    return;
  }
  WaveParameters *pars = (WaveParameters*) params;

  if(!reader.hasSection(string("Wave"))){
    return;
  }

  if(reader.hasParameter(string("Wave"),string("InitialConditions"))){
    string result = reader.readAsString(string("Wave"),string("InitialConditions"));
    if(result.compare("GAUSSIAN") == 0){
      pars->setInitialConditions(WaveParameters::GAUSSIAN);
    }
    else if(result.compare("FLAT") == 0){
      pars->setInitialConditions(WaveParameters::FLAT);
    }
    else {
      cout << "The value for parameter InitialConditions is out of range.\n";
    }
  }

  if(reader.hasParameter(string("Wave"),string("KOSigma"))){
    double result = reader.readAsDouble(string("Wave"),string("KOSigma"));
    if(result <= 1.000000e+00 && result >= 0.000000e+00){
      pars->setKOSigma(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

}

