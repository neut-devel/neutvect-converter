#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

#include "nvconv.h"

#include "NuHepMC/AttributeUtils.hxx"
#include "NuHepMC/make_writer.hxx"

#include <iostream>

std::vector<std::string> files_to_read;
std::string file_to_write;

std::string flux_file;
std::string flux_hist;

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
      << "\t-f <flux_file,flux_hist>     : ROOT flux histogram to use to\n"
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
      std::cout << "[INFO]: Assuming input flux histogram is in GeV."
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
        flux_hist = arg.substr(arg.find_first_of(',') + 1);
        std::cout << "[INFO]: Reading flux information from " << flux_file
                  << ":" << flux_hist << std::endl;
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

double GetFATX(TFile *fin, TChain &chin, NeutVect *nv, TH1 *&flux_histo,
               bool &isMonoE, int &beam_pid, double &flux_energy_to_MeV) {

  // reset it so that the caller knows if it comes back set, it was set by this
  // call
  flux_histo = nullptr;

  Long64_t ents = chin.GetEntries();
  chin.GetEntry(0);

  // check if it is mono-energetic
  double first_E = 0;
  isMonoE = true;
  for (Long64_t i = 0; i < std::min(1000ll, ents); ++i) {
    chin.GetEntry(i);
    if (!i) {
      first_E = nv->PartInfo(0)->fP.E();
      beam_pid = nv->PartInfo(0)->fPID;
    } else if (std::fabs(first_E - nv->PartInfo(0)->fP.E()) > 1E-6) {
      isMonoE = false;
    }
  }

  if (isMonoE) {
    double fatx = nv->Totcrs * 1E-2;
    std::cout << "[INFO]: Calculated FATX/event from Totcrs for assumed "
                 "mono-energetic file: "
              << fatx << " pb/Nucleon" << std::endl;
    return fatx;
  } else {
    std::cout << "[INFO]: Not Mono E, so cannot infer FATX from first event."
              << std::endl;
  }

  // if we have a flux file then we can build it
  if (flux_file.length() && flux_hist.length()) {

    std::unique_ptr<TFile> flux_fin(TFile::Open(flux_file.c_str(), "READ"));

    if (!flux_fin || !flux_fin->IsOpen() || flux_fin->IsZombie()) {
      std::cout << "[ERROR]: When trying to determine flux-average total cross "
                   "section, failed to open "
                << flux_file << std::endl;
      abort();
    }

    flux_histo = flux_fin->Get<TH1>(flux_hist.c_str());
    if (!flux_histo) {
      std::cout << "[ERROR]: When trying to determine flux-average total cross "
                   "section, failed to get "
                << flux_hist << " from " << flux_file << std::endl;
      abort();
    }

    flux_histo = static_cast<TH1 *>(flux_histo->Clone("fluxhisto_clone"));
    flux_histo->SetDirectory(nullptr);

    std::unique_ptr<TH1> xsechisto(
        static_cast<TH1 *>(flux_histo->Clone("xsechisto")));
    xsechisto->SetDirectory(nullptr);
    xsechisto->Reset();
    std::unique_ptr<TH1> entryhisto(
        static_cast<TH1 *>(flux_histo->Clone("entryhisto")));
    entryhisto->SetDirectory(nullptr);
    entryhisto->Reset();

    for (Long64_t i = 0; i < ents; ++i) {
      chin.GetEntry(i);
      double E_flux_units = nv->PartInfo(0)->fP.E() * (flux_in_GeV ? 1E-3 : 1);
      xsechisto->Fill(E_flux_units, nv->Totcrs);
      entryhisto->Fill(E_flux_units);
    }

    if (flux_in_GeV) {
      flux_energy_to_MeV = 1E3;
    } else {
      flux_energy_to_MeV = 1;
    }

    xsechisto->Divide(entryhisto.get());

    std::unique_ptr<TH1> ratehisto(
        static_cast<TH1 *>(xsechisto->Clone("ratehisto")));
    ratehisto->SetDirectory(nullptr);

    ratehisto->Multiply(flux_histo);

    double fatx = 1E-2 * (ratehisto->Integral() / flux_histo->Integral());

    flux_fin->Close();

    std::cout << "[INFO]: Calculated FATX/event from input file file as: "
              << fatx << " pb/Nucleon" << std::endl;
    return fatx;
  } else {
    std::cout << "[INFO]: Was not passed a flux file to calculat FATX"
              << std::endl;
  }

  // Have the histos already
  if (fin->Get<TH1>("ratehisto") && fin->Get<TH1>("fluxhisto")) {
    TH1 *ratehisto = fin->Get<TH1>("ratehisto");
    flux_histo = static_cast<TH1 *>(
        fin->Get<TH1>("fluxhisto")->Clone("fluxhisto_clone"));
    flux_histo->SetDirectory(nullptr);

    double fatx = 1E-2 * (ratehisto->Integral() / flux_histo->Integral());

    std::cout
        << "[INFO]: Calculated FATX/event from histograms in input file as: "
        << fatx << " pb/Nucleon" << std::endl;
    return fatx;
  } else {
    std::cout << "[INFO]: Cannot find ratehisto and fluxhisto in input file."
              << std::endl;
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
  auto branch_status = chin.SetBranchAddress("vectorbranch", &nv);

  if (skip >= ents) {
    std::cout << "Skipping " << skip << ", but only have " << ents
              << " in the input file." << std::endl;
    return 1;
  }

  Long64_t ents_to_run = std::min(ents, skip + nmaxevents);
  Long64_t ents_to_process = ents_to_run - skip;

  auto first_file = std::unique_ptr<TFile>(
      TFile::Open(files_to_read.front().c_str(), "READ"));

  TH1 *flux_histo = nullptr;
  bool isMonoE = false;
  int beam_pid = 0;

  double flux_energy_to_MeV = 1E3;
  double fatx = GetFATX(first_file.get(), chin, nv, flux_histo, isMonoE,
                        beam_pid, flux_energy_to_MeV);
  first_file->Close();
  first_file = nullptr;

  chin.GetEntry(0);

  int molecule_A = nv->TargetA;
  int molecule_H = nv->TargetH;

  auto gri = BuildRunInfo(ents_to_run, fatx, flux_histo, isMonoE, beam_pid,
                          flux_energy_to_MeV);

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

    auto hepev = ToGenEvent(nv, gri);

    hepev.set_event_number(i);
    NuHepMC::add_attribute(hepev, "ifile.name", fname);
    NuHepMC::add_attribute(hepev, "ifile.entry", fentry++);

    output->write_event(hepev);
  }
  std::cout << "\rConverting " << ents_to_process << "/" << ents_to_process
            << std::endl;

  output->close();
}