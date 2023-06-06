#include "TFile.h"
#include "TTree.h"

#include "neutvect.h"

#include "nvconv.h"

#include "HepMC3/WriterAscii.h"

#include <iostream>

int main(int argc, char const *argv[]) {

  if(argc <= 2){
    std::cout << "[ERROR]: Expected at least two arguments." << std::endl;
    return 1;
  }

  TFile *fin = TFile::Open(argv[1], "READ");

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

  HepMC3::WriterAscii output(argv[2]);
  if (output.failed()) {
    return 2;
  }

  for (Long64_t i = 0; i < ents; ++i) {
    tin->GetEntry(i);
    output.write_event(ToGenEvent(nv));
  }

  output.close();
}