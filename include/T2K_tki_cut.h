#pragma once
#include "tkigeneral.h"

#include <TLorentzVector.h>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

class TKIEvent;

TLorentzVector get_full_hadron_t2k(const TKIEvent &e);

ROOT::RDF::RNode DoTKICut_T2K(ROOT::RDF::RNode df);