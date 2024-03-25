#pragma once
#include "tkigeneral.h"

#include <TLorentzVector.h>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

class TKIEvent;

TLorentzVector get_full_hadron_MINERvAPi0(const TKIEvent &e);

ROOT::RDF::RNode DoTKICut_MINERvA(ROOT::RDF::RNode df) ;