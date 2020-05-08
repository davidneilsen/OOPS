#include <maxwellparser.h>
#include <iostream>

using namespace std;

void MaxwellParser::updateParameters(string fname, Parameters* params){
  if(!checkId(params)){
    return;
  }
  reader.clearData();
  ParamReader::ParamResult result = reader.readFile(fname);
  if(result != ParamReader::SUCCESS){
    cout << "An error occurred while trying to read " << fname << ".\n";
    return;
  }
  MaxwellParameters *pars = (MaxwellParameters*) params;

  if(!reader.hasSection(string("Maxwell"))){
    return;
  }

  if(reader.hasParameter(string("Maxwell"),string("npoints"))){
    int result = reader.readAsInt(string("Maxwell"),string("npoints"));
    if(result <= 50001 && result >= 5){
      pars->setnpoints(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("xmin"))){
    double result = reader.readAsDouble(string("Maxwell"),string("xmin"));
    if(result <= 1.000000e+03 && result >= -1.000000e+03){
      pars->setxmin(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("xmax"))){
    double result = reader.readAsDouble(string("Maxwell"),string("xmax"));
    if(result <= 1.000000e+03 && result >= -1.000000e+03){
      pars->setxmax(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("tmax"))){
    double result = reader.readAsDouble(string("Maxwell"),string("tmax"));
    if(result <= 1.000000e+05 && result >= 0.000000e+00){
      pars->settmax(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("cfl"))){
    double result = reader.readAsDouble(string("Maxwell"),string("cfl"));
    if(result <= 1.000000e+00 && result >= 0.000000e+00){
      pars->setcfl(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("id_gauss_center"))){
    double result = reader.readAsDouble(string("Maxwell"),string("id_gauss_center"));
    if(result <= 1.000000e+03 && result >= -1.000000e+03){
      pars->setid_gauss_center(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("id_gauss_width"))){
    double result = reader.readAsDouble(string("Maxwell"),string("id_gauss_width"));
    if(result <= 1.000000e+01 && result >= 1.000000e-03){
      pars->setid_gauss_width(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

  if(reader.hasParameter(string("Maxwell"),string("KOSigma"))){
    double result = reader.readAsDouble(string("Maxwell"),string("KOSigma"));
    if(result <= 1.000000e+00 && result >= 0.000000e+00){
      pars->setKOSigma(result);
    }
    else{
      cout << "Parameter %s out of range.\n";
    }
  }

}

