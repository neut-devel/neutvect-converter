#pragma once

#include "neutvect.h"

#include "HepMC3/GenEvent.h"

#include "TH1.h"

namespace NuHepMC {

namespace ParticleStatus {
// NEUT Extended codes
namespace NEUT {
const int PauliBlocked = 31;
const int UnderwentFSI = 41;
const int SecondaryInteraction = 99;
} // namespace NEUT
} // namespace ParticleStatus

} // namespace NuHepMC

std::shared_ptr<HepMC3::GenRunInfo>
BuildRunInfo(int nevents, double flux_averaged_total_cross_section,
             TH1 *&flux_histo, bool &isMonoE, int beam_pid,
             double flux_to_MeV = 1);
HepMC3::GenEvent ToGenEvent(NeutVect *nv,
                            std::shared_ptr<HepMC3::GenRunInfo> gri);