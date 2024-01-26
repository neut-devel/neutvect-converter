#pragma once

#include "neutvect.h"

#include "TChain.h"
#include "TH1.h"

#include <memory>
#include <optional>
#include <utility>

namespace nvconv {

bool isMono(TChain &chin, NeutVect *nv, Long64_t ntocheck = 1000ll);

std::unique_ptr<TH1> GetHistFromFile(std::string const &flux_file,
                                     std::string const &flux_histname);
std::optional<double> GetFATXFromFluxHist(TChain &chin, NeutVect *nv,
                                          std::unique_ptr<TH1> const &flux_hist,
                                          bool flux_in_GeV);

std::pair<std::unique_ptr<TH1>, std::unique_ptr<TH1>>
GetFluxRateHistPairFromChain(TChain &chin);

} // namespace nvconv