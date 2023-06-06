#pragma once

#include "neutvect.h"
#include "HepMC3/GenEvent.h"

namespace nvconverter {

// HepMC3 Status Codes
static int const kISStatus = 4;
static int const kFSStatus = 1;

static int const kTgtNucleusStatus = 11;
static int const kOther = 12;

// 1: Decayed to the other particle
static int const kFSIDecayed = 13;
// 3: Absorped
static int const kFSIAbsorpbed = 14;
// 4: Charge exchanged
static int const kFSIChargeExchanged = 15;
// 5: Pauli blocked
static int const kPauliBlocked = 16;

} // namespace nvconverter

HepMC3::GenEvent ToGenEvent(NeutVect *nv,
                            double flux_averaged_total_cross_section = 1);