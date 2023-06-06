#pragma once

#include "HepMC3/GenEvent.h"
#include "neutvect.h"

namespace NEUTNuHepMC {

namespace VertexStatus {
// NuHepMC standard vertex codes
const int kPrimaryVertex = 1;
const int kNuclearTargetVertex = 2;

// NEUT Extended codes
const int kFSIVertex = 3;
} // namespace VertexStatus

namespace ParticleStatus {
// HepMC3 standard vertex codes
const int kUndecayedPhysicalParticle = 1;
const int kDecayedParticle = 2;
const int kDocumentationLine = 3;
const int kIncomingBeamParticle = 4;

// NuHepMC standard vertex codes
const int kTargetParticle = 11;

// P.C.1
const int kStruckNucleon = 21;

// NEUT Extended codes
const int kPauliBlocked = 31;
const int kUnderwentFSI = 41;
const int kSecondaryInteraction = 99;
} // namespace ParticleStatus

} // namespace NEUTNuHepMC

std::shared_ptr<HepMC3::GenRunInfo>
BuildRunInfo(int nevents, double flux_averaged_total_cross_section = 1);
HepMC3::GenEvent ToGenEvent(NeutVect *nv, std::shared_ptr<HepMC3::GenRunInfo> gri);