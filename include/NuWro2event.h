#pragma once

#include "event1.h"
#include "tkievent.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>

NeutrinoEvent NuWro2event(event &e1);
NeutrinoEvent NuWro2event_nofsi(event &e1);

ROOT::RDF::RNode NuWroPrepare(ROOT::RDF::RNode df, bool fsi = true);

