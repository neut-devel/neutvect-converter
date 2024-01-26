#include "nvfatxtools.h"

#include "TFile.h"

namespace nvconv {

bool isMono(TChain &chin, NeutVect *nv, Long64_t ntocheck) {
  Long64_t ents = chin.GetEntries();
  chin.GetEntry(0);

  double first_E = 0;
  for (Long64_t i = 0; i < std::min(ntocheck, ents); ++i) {
    chin.GetEntry(i);
    if (!i) {
      first_E = nv->PartInfo(0)->fP.E();
    } else if (std::fabs(first_E - nv->PartInfo(0)->fP.E()) > 1E-6) {
      return false;
    }
  }
  return true;
}

std::unique_ptr<TH1> GetHistFromFile(std::string const &flux_file,
                                     std::string const &flux_histname) {

  if(!flux_file.length() || !flux_histname.length()){
    return nullptr;
  }

  std::unique_ptr<TFile> flux_fin(TFile::Open(flux_file.c_str(), "READ"));

  if (!flux_fin || !flux_fin->IsOpen() || flux_fin->IsZombie()) {
    std::cout << "[ERROR]: When trying to determine flux-average total cross "
                 "section, failed to open "
              << flux_file << std::endl;
    return nullptr;
  }

  TH1 *flux_hist = flux_fin->Get<TH1>(flux_histname.c_str());
  if (!flux_hist) {
    std::cout << "[ERROR]: When trying to determine flux-average total cross "
                 "section, failed to get "
              << flux_hist << " from " << flux_file << std::endl;
    return nullptr;
  }

  std::unique_ptr<TH1> fhc(
      static_cast<TH1 *>(flux_hist->Clone("fluxhisto_clone")));
  fhc->SetDirectory(nullptr);

  flux_fin->Close();

  return fhc;
}

std::optional<double> GetFATXFromFluxHist(TChain &chin, NeutVect *nv,
                                          std::unique_ptr<TH1> const &flux_hist,
                                          bool flux_in_GeV) {

  if (!flux_hist) {
    return std::optional<double>();
  }

  std::unique_ptr<TH1> xsechisto(
      static_cast<TH1 *>(flux_hist->Clone("xsechisto")));
  xsechisto->SetDirectory(nullptr);
  xsechisto->Reset();
  std::unique_ptr<TH1> entryhisto(
      static_cast<TH1 *>(flux_hist->Clone("entryhisto")));
  entryhisto->SetDirectory(nullptr);
  entryhisto->Reset();

  Long64_t ents = chin.GetEntries();
  for (Long64_t i = 0; i < ents; ++i) {
    chin.GetEntry(i);
    double E_flux_units = nv->PartInfo(0)->fP.E() * (flux_in_GeV ? 1E-3 : 1);
    xsechisto->Fill(E_flux_units, nv->Totcrs);
    entryhisto->Fill(E_flux_units);
  }

  xsechisto->Divide(entryhisto.get());

  std::unique_ptr<TH1> ratehisto(
      static_cast<TH1 *>(xsechisto->Clone("ratehisto")));
  ratehisto->SetDirectory(nullptr);

  ratehisto->Multiply(flux_hist.get());

  double fatx = 1E-2 * (ratehisto->Integral() / flux_hist->Integral());

  return fatx;
}

std::pair<std::unique_ptr<TH1>, std::unique_ptr<TH1>>
GetFluxRateHistPairFromChain(TChain &chin) {
  chin.GetEntry(0);

  if (chin.GetFile()->Get<TH1>("ratehisto") &&
      chin.GetFile()->Get<TH1>("fluxhisto")) {
    std::unique_ptr<TH1> rhc(static_cast<TH1 *>(
        chin.GetFile()->Get<TH1>("ratehisto")->Clone("ratehisto_clone")));
    rhc->SetDirectory(nullptr);
    std::unique_ptr<TH1> fhc(static_cast<TH1 *>(
        chin.GetFile()->Get<TH1>("fluxhisto")->Clone("fluxhisto_clone")));
    fhc->SetDirectory(nullptr);
    return {std::move(rhc), std::move(fhc)};
  }
  return {nullptr, nullptr};
}

} // namespace nvconv