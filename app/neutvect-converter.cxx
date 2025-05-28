#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

#include "nvconv.h"
#include "nvfatxtools.h"

#include "NuHepMC/AttributeUtils.hxx"
#include "NuHepMC/make_writer.hxx"

#include <iostream>

std::vector<std::string> files_to_read;
std::string file_to_write;

std::string flux_file = "";
std::string flux_histname = "";

bool flux_in_GeV = true;
double monoE = 0;
Long64_t skip = 0;

Long64_t nmaxevents = std::numeric_limits<Long64_t>::max();

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0] << "\n"
      << "\t-i <nv.root> [nv2.root ...]  : neutvect file to read\n"
      << "\t-N <NMax>                    : Process at most <NMax> events\n"
      << "\t-o <neut.hepmc3>             : hepmc3 file to write\n"
      << "\t-f <flux_file,flux_histname>     : ROOT flux histogram to use to\n"
      << "\t-M                           : -f argument should be interpreted "
         "as being in MeV\n"
      << "\t-s <N>                       : Skip <N>." << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?" || std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-M") {
      flux_in_GeV = false;
      std::cout << "[INFO]: Assuming input flux histogram is in MeV."
                << std::endl;
    } else if ((opt + 1) < argc) {
      if (std::string(argv[opt]) == "-i") {
        while (((opt + 1) < argc) && (argv[opt + 1][0] != '-')) {
          files_to_read.push_back(argv[++opt]);
          std::cout << "[INFO]: Reading from " << files_to_read.back()
                    << std::endl;
        }
      } else if (std::string(argv[opt]) == "-N") {
        nmaxevents = std::stol(argv[++opt]);
        std::cout << "[INFO]: Processing at most " << nmaxevents << " events."
                  << std::endl;
      } else if (std::string(argv[opt]) == "-s") {
        skip = std::stol(argv[++opt]);
        std::cout << "[INFO]: Skipping " << skip << " events before processing."
                  << std::endl;
      } else if (std::string(argv[opt]) == "-o") {
        file_to_write = argv[++opt];
      } else if (std::string(argv[opt]) == "-f") {
        std::string arg = argv[++opt];
        flux_file = arg.substr(0, arg.find_first_of(','));
        flux_histname = arg.substr(arg.find_first_of(',') + 1);
        std::cout << "[INFO]: Reading flux information from " << flux_file
                  << ":" << flux_histname << std::endl;
      } else {
        std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

double GetFATX(TChain &chin, NeutVect *nv, std::unique_ptr<TH1> &flux_hist,
               bool &isMonoE, int &beam_pid, double &flux_energy_to_MeV) {

  // reset it so that the caller knows if it comes back set, it was set by this
  // call
  flux_hist = nullptr;

  chin.GetEntries();
  chin.GetEntry(0);
  beam_pid = nv->PartInfo(0)->fPID;

  isMonoE = nvconv::isMono(chin, nv);
  if (isMonoE) {
    chin.GetEntry(0);
    double fatx = nv->Totcrs * 1E-2;
    std::cout << "[INFO]: Calculated FATX from Totcrs for assumed "
                 "mono-energetic file: "
              << fatx << " pb/Nucleon" << std::endl;
    return fatx;
  } else {
    std::cout
        << "[INFO]: Not mono-energetic, so cannot infer FATX from first event."
        << std::endl;
  }

  flux_hist = nvconv::GetHistFromFile(flux_file, flux_histname);

  // if we have a flux file then we can build it
  if (flux_hist) {
    auto fatx_opt =
        nvconv::GetFATXFromFluxHist(chin, nv, flux_hist, flux_in_GeV);
    if (fatx_opt) {
      if (flux_in_GeV) {
        flux_energy_to_MeV = 1E3;
      } else {
        flux_energy_to_MeV = 1;
      }
      std::cout << "[INFO]: Calculated FATX from input file file as: "
                << fatx_opt.value() << " pb/Nucleon" << std::endl;
      return fatx_opt.value();
    }
  }

  auto frpair = nvconv::GetFluxRateHistPairFromChain(chin);
  if (frpair.second) {

    double fatx = 1E-2 * (frpair.first->Integral() / frpair.second->Integral());

    flux_hist = std::move(frpair.second);

    std::cout
        << "[INFO]: Calculated FATX from histograms in input file as: 1E-2 * "
        << frpair.first->Integral() << "/" << flux_hist->Integral() << " = "
        << fatx << " pb/Nucleon" << std::endl;

    return fatx;
  }

  return 1;
}

int main(int argc, char const *argv[]) {

  handleOpts(argc, argv);

  if (!files_to_read.size() || !file_to_write.length()) {
    std::cout << "[ERROR]: Expected -i and -o arguments." << std::endl;
    return 1;
  }

  TChain chin("neuttree");

  for (auto const &ftr : files_to_read) {
    if (!chin.Add(ftr.c_str(), 0)) {
      std::cout << "[ERROR]: Failed to find tree: \"neuttree\" in file: \""
                << ftr << "\"." << std::endl;
      return 1;
    }
  }

  Long64_t ents = chin.GetEntries();
  chin.GetEntry(
      0); // need to do this before opening the other file or... kablamo

  NeutVect *nv = nullptr;
  chin.SetBranchAddress("vectorbranch", &nv);

  if (skip >= ents) {
    std::cout << "Skipping " << skip << ", but only have " << ents
              << " in the input file." << std::endl;
    return 1;
  }

  Long64_t ents_to_run = std::min(ents, skip + nmaxevents);
  Long64_t ents_to_process = ents_to_run - skip;

  auto first_file = std::unique_ptr<TFile>(
      TFile::Open(files_to_read.front().c_str(), "READ"));

  std::unique_ptr<TH1> flux_histo = nullptr;
  bool isMonoE = false;
  int beam_pid = 0;

  double flux_energy_to_MeV = 1E3;
  double fatx =
      GetFATX(chin, nv, flux_histo, isMonoE, beam_pid, flux_energy_to_MeV);
  first_file->Close();
  first_file = nullptr;

  chin.GetEntry(0);

  int molecule_A = nv->TargetA;
  int molecule_H = nv->TargetH;

  auto gri = nvconv::BuildRunInfo(ents_to_run, fatx, flux_histo, isMonoE,
                                  beam_pid, flux_energy_to_MeV);

  std::unique_ptr<HepMC3::Writer> output(
      NuHepMC::Writer::make_writer(file_to_write, gri));

  if (output->failed()) {
    return 2;
  }

  Long64_t fentry = 0;
  TUUID fuid = chin.GetFile()->GetUUID();
  std::string fname = chin.GetFile()->GetName();
  for (Long64_t i = 0; i < ents_to_run; ++i) {
    if (i >= ents) {
      break;
    }
    chin.GetEntry(i);

    if (i == 0) {
      fuid = chin.GetFile()->GetUUID();
      fname = chin.GetFile()->GetName();
    } else if (chin.GetFile()->GetUUID() != fuid) {
      fuid = chin.GetFile()->GetUUID();
      fname = chin.GetFile()->GetName();
      fentry = 0;
    }

    if (i < skip) { // we have to manually skip them to work out the correct
                    // file indexing
      fentry++;
      continue;
    }

    if ((molecule_A != nv->TargetA) || (molecule_H != nv->TargetH)) {
      std::cout << "neutvect-converter cannot currently convert to NuHepMC for "
                   "multi-target event vectors."
                << std::endl;
      return 1;
    }

    if (i && (ents_to_process / 100) && !(i % (ents_to_process / 100))) {
      std::cout << "\rConverting " << i << "/" << ents_to_process << std::flush;
    }

    auto hepev = nvconv::ToGenEvent(nv, gri);

    hepev->set_event_number(i);
    NuHepMC::add_attribute(*hepev, "ifile.name", fname);
    NuHepMC::add_attribute(*hepev, "ifile.entry", fentry++);

    output->write_event(*hepev);
  }
  std::cout << "\rConverting " << ents_to_process << "/" << ents_to_process
            << std::endl;

  output->close();
}