#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "neutvect.h"

#include "nvconv.h"

#include "HepMC3/WriterAscii.h"

#include <iostream>

std::string file_to_read;
std::string file_to_write;

std::string flux_file;
std::string flux_hist;

bool flux_in_GeV = false;
double monoE = 0;

Long64_t nmaxevents = std::numeric_limits<Long64_t>::max();

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0] << "\n"
      << "\t-i <neutvect.root>       : neutvect file to read\n"
      << "\t-N <NMax>                : Process at most <NMax> events\n"
      << "\t-o <neut.hepmc3>         : hepmc3 file to write\n"
      << "\t-f <flux_file,flux_hist> : ROOT flux histogram to use to\n"
      << "\t-G                       : -f argument should be interpreted as "
         "being in GeV\n"
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
    } else if ((opt + 1) < argc) {
      if (std::string(argv[opt]) == "-i") {
        file_to_read = argv[++opt];
        std::cout << "[INFO]: Reading from " << file_to_read << std::endl;
      } else if (std::string(argv[opt]) == "-N") {
        nmaxevents = std::stol(argv[++opt]);
        std::cout << "[INFO]: Processing at most " << nmaxevents << " events."
                  << std::endl;
      } else if (std::string(argv[opt]) == "-o") {
        file_to_write = argv[++opt];
        std::cout << "[INFO]: Writing to " << file_to_write << std::endl;
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

double GetFATX(TFile *fin, TTree *tin, NeutVect *nv) {

  Long64_t ents = tin->GetEntries();

  // check if it is mono-energetic
  double first_E = 0;
  bool isMonoE = true;
  for (Long64_t i = 0; i < std::min(1000ll, ents); ++i) {
    tin->GetEntry(i);
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
      tin->GetEntry(i);
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

  if (!file_to_read.length() || !file_to_write.length()) {
    std::cout << "[ERROR]: Expected -i and -o arguments." << std::endl;
    return 1;
  }

  TFile *fin = TFile::Open(file_to_read.c_str(), "READ");
  if (!fin) {
    std::cout << "[ERROR]: Failed to read input file: " << argv[1] << std::endl;
    return 1;
  }

  TTree *tin = fin->Get<TTree>("neuttree");
  if (!tin) {
    std::cout << "[ERROR]: Failed to read input tree: neuttree, from file: "
              << argv[1] << std::endl;
    return 1;
  }

  NeutVect *nv = nullptr;
  tin->SetBranchAddress("vectorbranch", &nv);

  Long64_t ents = tin->GetEntries();
  Long64_t ents_to_run = std::min(ents, nmaxevents);

  double fatx = GetFATX(fin, tin, nv) * (double(ents_to_run) / double(ents));

  auto gri = BuildRunInfo(ents_to_run, fatx);
  HepMC3::WriterAscii output(file_to_write.c_str(), gri);

  tin->GetEntry(0);

  if (output.failed()) {
    return 2;
  }

  for (Long64_t i = 0; i < ents_to_run; ++i) {
    tin->GetEntry(i);
    if (i && (ents_to_run / 100) && !(i % (ents_to_run / 100))) {
      std::cout << "\rConverting " << i << "/" << ents_to_run << std::flush;
    }
    output.write_event(ToGenEvent(nv, gri));
  }
  std::cout << "\rConverting " << ents_to_run << "/" << ents_to_run
            << std::endl;

  output.close();
}