
#include "NuHepMC/CrossSectionUtils.hxx"

#include "TFile.h"
#include "TGraph.h"

int main(int argc, char const *argv[]) {
  auto splines = NuHepMC::CrossSection::BuildSplines(argv[1]);

  TFile out(argv[2], "RECREATE");

  for (auto const &beam_map : splines) {
    for (auto const &tgt_map : beam_map.second) {
      for (auto const &proc_splines : tgt_map.second) {

        TGraph g;
        size_t npoints = 0;
        for (auto sp : proc_splines.second.points) {
          g.SetPoint(npoints++, sp.first, sp.second);
        }

        std::string spl_name = "spline_nupid_";
        spl_name += std::to_string(beam_map.first);
        spl_name += "tgtpid_";
        spl_name += std::to_string(tgt_map.first);
        spl_name +=
            proc_splines.first == 0
                ? std::string("totxs")
                : (std::string("procid_") + std::to_string(proc_splines.first));

        out.WriteTObject(&g, spl_name.c_str());
      }
    }
  }

  out.Write();
  out.Close();
}