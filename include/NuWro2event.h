#pragma once

#include "event1.h"
#include "tkievent.h"

#include <TLorentzVector.h>
#include "ROOT/RDF/InterfaceUtils.hxx"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>


TKIEvent NuWro2event(event &e1);

ROOT::RDF::RNode NuWroPrepare(ROOT::RDF::RNode df);

ROOT::RDF::RNode CommonVariableDefine(ROOT::RDF::RNode df);
