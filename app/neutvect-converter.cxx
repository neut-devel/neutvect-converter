#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

#include "nvconv.h"

#include "HepMC3/WriterAscii.h"
#ifdef HEPMC3_USE_COMPRESSION
#include "HepMC3/WriterGZ.h"
#endif

#include <iostream>

std::vector<std::string> files_to_read;
std::string file_to_write;

std::string flux_file;
std::string flux_hist;
#ifdef HEPMC3_USE_COMPRESSION
bool WriteGZ = false;
#endif

bool flux_in_GeV = false;
double monoE = 0;

Long64_t nmaxevents = std::numeric_limits<Long64_t>::max();

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0] << "\n"
      << "\t-i <nv.root> [nv2.root ...]  : neutvect file to read\n"
      << "\t-N <NMax>                    : Process at most <NMax> events\n"
      << "\t-o <neut.hepmc3>             : hepmc3 file to write\n"
#ifdef HEPMC3_USE_COMPRESSION
      << "\t-z                           : write out in compressed ASCII\n"
#endif
      << "\t-f <flux_file,flux_hist>     : ROOT flux histogram to use to\n"
      << "\t-G                           : -f argument should be interpreted "
         "as being in GeV\n"
      << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?" || std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-G") {
      flux_in_GeV = true;
      std::cout << "[INFO]: Assuming input flux histogram is in GeV."
                << std::endl;
    } else if (std::string(argv[opt]) == "-z") {
      WriteGZ = true;
      std::cout << "[INFO]: Writing output compressed output file."
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
      } else if (std::string(argv[opt]) == "-o") {
        file_to_write = argv[++opt];
      } else if (std::string(argv[opt]) == "-f") {
        std::string arg = argv[++opt];
        flux_file = arg.substr(0, arg.find_first_of(','));
        flux_hist = arg.substr(arg.find_first_of(',') + 1);
        std::cout << "[INFO]: Reading flux information from " << flux_file
                  << ":" << flux_hist << std::endl;
      }
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

double GetFATX(TFile *fin, TChain &chin, NeutVect *nv) {

  Long64_t ents = chin.GetEntries();
  chin.GetEntry(0);

  // check if it is mono-energetic
  double first_E = 0;
  bool isMonoE = true;
  for (Long64_t i = 0; i < std::min(1000ll, ents); ++i) {
    chin.GetEntry(i);
    if (!i) {
      first_E = nv->PartInfo(0)->fP.E();
    } else if (std::fabs(first_E - nv->PartInfo(0)->fP.E()) > 1E-6) {
      isMonoE = false;
    }
  }

  if (isMonoE) {
    double fatx = nv->Totcrs * 1E-38 / double(ents);
    std::cout << "[INFO]: Calculated FATX/event from Totcrs for assumed "
                 "mono-energetic file: "
              << fatx << std::endl;
    return fatx;
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

    TH1 *fluxhisto = flux_fin->Get<TH1>(flux_hist.c_str());
    if (!fluxhisto) {
      std::cout << "[ERROR]: When trying to determine flux-average total cross "
                   "section, failed to get "
                << flux_hist << " from " << flux_file << std::endl;
      abort();
    }

    std::unique_ptr<TH1> xsechisto(
        static_cast<TH1 *>(fluxhisto->Clone("xsechisto")));
    xsechisto->SetDirectory(nullptr);
    xsechisto->Reset();
    std::unique_ptr<TH1> entryhisto(
        static_cast<TH1 *>(fluxhisto->Clone("entryhisto")));
    entryhisto->SetDirectory(nullptr);
    entryhisto->Reset();

    for (Long64_t i = 0; i < ents; ++i) {
      chin.GetEntry(i);
      double E = nv->PartInfo(0)->fP.E() * (flux_in_GeV ? 1E-3 : 1);
      xsechisto->Fill(E, nv->Totcrs);
      entryhisto->Fill(E);
    }

    xsechisto->Divide(entryhisto.get());

    std::unique_ptr<TH1> ratehisto(
        static_cast<TH1 *>(xsechisto->Clone("ratehisto")));
    ratehisto->SetDirectory(nullptr);

    ratehisto->Multiply(fluxhisto);

    double fatx = 1E-38 * (ratehisto->Integral() /
                           (fluxhisto->Integral() * double(ents)));

    flux_fin->Close();

    std::cout << "[INFO]: Calculated FATX/event from input file file as: "
              << fatx << std::endl;
    return fatx;
  }

  // Have the histos already
  if (fin->Get<TH1D>("ratehisto") && fin->Get<TH1D>("fluxhisto")) {
    TH1D *ratehisto = fin->Get<TH1D>("ratehisto");
    TH1D *fluxhisto = fin->Get<TH1D>("fluxhisto");

    double fatx = 1E-38 * (ratehisto->Integral() /
                           (fluxhisto->Integral() * double(ents)));

    std::cout
        << "[INFO]: Calculated FATX/event from histograms in input file as: "
        << fatx << std::endl;
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

#ifdef HEPMC3_USE_COMPRESSION
  if (WriteGZ && (file_to_write.substr(file_to_write.size() - 4, 3) !=
                  std::string(".gz"))) {
    file_to_write = file_to_write + ".gz";
  }
  std::cout << "[INFO]: Writing to " << file_to_write << std::endl;
#endif

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

  Long64_t ents_to_run = std::min(ents, nmaxevents);

  auto first_file = std::unique_ptr<TFile>(
      TFile::Open(files_to_read.front().c_str(), "READ"));

  double fatx = GetFATX(first_file.get(), chin, nv) *
                (double(ents)/double(ents_to_run));
  first_file->Close();
  first_file = nullptr;

  auto gri = BuildRunInfo(ents_to_run, fatx);
  HepMC3::Writer *output =
#ifdef HEPMC3_USE_COMPRESSION
      WriteGZ ? static_cast<HepMC3::Writer *>(
                    new HepMC3::WriterGZ<HepMC3::WriterAscii>(
                        file_to_write.c_str(), gri))
              :
#endif
              static_cast<HepMC3::Writer *>(
                  new HepMC3::WriterAscii(file_to_write.c_str(), gri));

  chin.GetEntry(0);

  if (output->failed()) {
    return 2;
  }

  for (Long64_t i = 0; i < ents_to_run; ++i) {
    chin.GetEntry(i);
    if (i && (ents_to_run / 100) && !(i % (ents_to_run / 100))) {
      std::cout << "\rConverting " << i << "/" << ents_to_run << std::flush;
    }
    output->write_event(ToGenEvent(nv, gri));
  }
  std::cout << "\rConverting " << ents_to_run << "/" << ents_to_run
            << std::endl;

  output->close();
}