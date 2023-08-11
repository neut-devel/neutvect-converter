#include "nvconv.h"

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include <utility>

using namespace NEUTNuHepMC;

// Lazy way of choosing the right attribute type via TMP
template <typename T> struct attr_traits {};

template <> struct attr_traits<int> {
  typedef HepMC3::IntAttribute type;
};

template <> struct attr_traits<std::vector<int>> {
  typedef HepMC3::VectorIntAttribute type;
};

template <> struct attr_traits<double> {
  typedef HepMC3::DoubleAttribute type;
};

template <> struct attr_traits<std::vector<double>> {
  typedef HepMC3::VectorDoubleAttribute type;
};

template <> struct attr_traits<std::string> {
  typedef HepMC3::StringAttribute type;
};

template <size_t N> struct attr_traits<char[N]> {
  typedef HepMC3::StringAttribute type;
};

template <> struct attr_traits<std::vector<std::string>> {
  typedef HepMC3::VectorStringAttribute type;
};

template <typename T>
void add_attribute(HepMC3::GenEvent &ge, std::string const &name,
                   T const &val) {
  ge.add_attribute(name, std::make_shared<typename attr_traits<T>::type>(val));
}

template <typename T>
void add_attribute(std::shared_ptr<HepMC3::GenRunInfo> gri,
                   std::string const &name, T const &val) {
  gri->add_attribute(name,
                     std::make_shared<typename attr_traits<T>::type>(val));
}

static const double GeV = 1E-3;
static const double MeV = 1;

std::shared_ptr<HepMC3::GenRunInfo>
BuildRunInfo(int nevents, double flux_averaged_total_cross_section) {

  // G.R.1 Valid GenRunInfo
  auto run_info = std::make_shared<HepMC3::GenRunInfo>();

  // G.R.2 NuHepMC Version
  add_attribute(run_info, "NuHepMC.Version.Major", 0);
  add_attribute(run_info, "NuHepMC.Version.Minor", 1);
  add_attribute(run_info, "NuHepMC.Version.Patch", 0);

  // G.R.3 Generator Identification
  run_info->tools().emplace_back(HepMC3::GenRunInfo::ToolInfo{
      "NEUT", NEUT_VERSION_STR,
      "https://doi.org/10.1140/epjs/s11734-021-00287-7"});
  run_info->tools().emplace_back(
      HepMC3::GenRunInfo::ToolInfo{"neutvect-converter", PROJECT_VERSION_STR,
                                   "github.com/neut-dev/neutvect-converter"});

  // G.R.4 Process Metadata
  // The reason for the second tuple entry is that this is copied wholesale from
  // NEUT6 which needs this map
  const std::vector<std::tuple<std::string, int, int>>
      ChannelNameIndexModeMapping{
          // from nemodsel.F on 2020/10/21
          // DATA MODNEU /  1, 11, 12, 13, 21, 31, 32, 33, 34, 41,
          // &              51, 51, 52, 16, 36, 22, 42, 43, 23, 44,
          // &              45,  0, 26, 46, 17, 38, 39, 2,  15, 35/

          {"CC_QE_nu", 1, 1},

          {"CC_RES_ppi+_nu", 2, 11},
          {"CC_RES_ppi0_nu", 3, 12},
          {"CC_RES_npi+_nu", 4, 13},

          {"CC_multi_pi_nu", 5, 21},

          {"NC_RES_npi0_nu", 6, 31},
          {"NC_RES_ppi0_nu", 7, 32},
          {"NC_RES_ppi-_nu", 8, 33},
          {"NC_RES_npi+_nu", 9, 34},

          {"NC_multi_pi_nu", 10, 41},

          {"NC_elastic_free_p_nu", 11, 51},
          {"NC_elastic_bound_p_nu", 12, 51},
          {"NC_elastic_n_nu", 13, 52},

          {"CC_COH_nu", 14, 16},

          {"NC_COH_nu", 15, 36},

          {"CC_eta_nu", 16, 22},

          {"NC_eta_n_nu", 17, 42},
          {"NC_eta_p_nu", 18, 43},

          {"CC_kaon_nu", 19, 23},

          {"NC_kaon_n_nu", 20, 44},
          {"NC_kaon_p_nu", 21, 45},

          {"CC_DIS_nu", 23, 26},

          {"NC_DIS_nu", 24, 46},

          {"CC_1gamma_nu", 25, 17},

          {"NC_1gamma_n_nu", 26, 38},
          {"NC_1gamma_p_nu", 27, 39},

          {"CC_2p2h_nu", 28, 2},

          {"CC_DIF_nu", 29, 15},

          {"NC_DIF_nu", 30, 35},

          // from nemodsel.F on 2020/10/21
          //  DATA MODNEUB/ -1,-11,-12,-13,-21,-31,-32,-33,-34,-41,
          // &              -1,-51,-51,-52,-16,-36,-22,-42,-43,-23,
          // &              -44,-45,-26,-46,-17,-38,-39,-2,-15,-35/

          {"CC_QE_free_proton_nubar", 1, -1},

          {"CC_RES_npi-_nubar", 2, -11},
          {"CC_RES_ppi0_nubar", 3, -12},
          {"CC_RES_ppi-_nubar", 4, -13},

          {"CC_multi_pi_nubar", 5, -21},

          {"NC_RES_npi0_nubar", 6, -31},
          {"NC_RES_ppi0_nubar", 7, -32},
          {"NC_RES_ppi-_nubar", 8, -33},
          {"NC_RES_npi+_nubar", 9, -34},

          {"NC_multi_pi_nubar", 10, -41},

          {"CC_QE_bound_proton_nubar", 11, -1},

          {"NC_elastic_free_p_nubar", 12, -51},
          {"NC_elastic_bound_p_nubar", 13, -51},
          {"NC_elastic_n_nubar", 14, -52},

          {"CC_COH_nubar", 15, -16},

          {"NC_COH_nubar", 16, -36},

          {"CC_eta_nubar", 17, -22},

          {"NC_eta_n_nubar", 18, -42},
          {"NC_eta_p_nubar", 19, -43},

          {"CC_kaon_nubar", 20, -23},

          {"NC_kaon_n_nubar", 21, -44},
          {"NC_kaon_p_nubar", 22, -45},

          {"CC_DIS_nubar", 23, -26},

          {"NC_DIS_nubar", 24, -46},

          {"CC_1gamma_nubar", 25, -17},

          {"NC_1gamma_n_nubar", 26, -38},
          {"NC_1gamma_p_nubar", 27, -39},

          {"CC_2p2h_nubar", 28, -2},

          {"CC_DIF_nubar", 29, -15},

          {"NC_DIF_nubar", 30, -35},

      };

  std::vector<int> pids;
  for (auto const &pid : ChannelNameIndexModeMapping) {
    pids.push_back(std::get<2>(pid));
    add_attribute(run_info,
                  "NuHepMC.ProcessInfo[" + std::to_string(std::get<2>(pid)) +
                      "].Name",
                  std::get<0>(pid));
    add_attribute(run_info,
                  "NuHepMC.ProcessInfo[" + std::to_string(std::get<2>(pid)) +
                      "].Description",
                  "");
  }

  add_attribute(run_info, "NuHepMC.ProcessIDs", pids);

  // G.R.5 Vertex Status Metadata
  add_attribute(run_info,
                "NuHepMC.VertexStatusInfo[" +
                    std::to_string(VertexStatus::kPrimaryVertex) + "].Name",
                "PrimaryVertex");
  add_attribute(run_info,
                "NuHepMC.VertexStatusInfo[" +
                    std::to_string(VertexStatus::kPrimaryVertex) +
                    "].Description",
                "The neutrino hard-scatter vertex");
  add_attribute(run_info,
                "NuHepMC.VertexStatusInfo[" +
                    std::to_string(VertexStatus::kFSIVertex) + "].Name",
                "FSIVertex");
  add_attribute(run_info,
                "NuHepMC.VertexStatusInfo[" +
                    std::to_string(VertexStatus::kFSIVertex) + "].Description",
                "A single vertex representing the cascade");

  add_attribute(
      run_info, "NuHepMC.VertexStatusIDs",
      std::vector<int>{VertexStatus::kPrimaryVertex, VertexStatus::kFSIVertex});

  // G.R.6 Particle Status Metadata
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kUndecayedPhysicalParticle) +
                    "].Name",
                "UndecayedPhysicalParticle");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kUndecayedPhysicalParticle) +
                    "].Description",
                "Physical final state particles produced by this simulation");

  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kDecayedParticle) + "].Name",
                "DecayedParticle");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kDecayedParticle) +
                    "].Description",
                "Particle was decayed by the simulation");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kDocumentationLine) +
                    "].Name",
                "DocumentationLine");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kDocumentationLine) +
                    "].Description",
                "Documentation line, not considered a real particle");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kIncomingBeamParticle) +
                    "].Name",
                "IncomingBeamParticle");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kIncomingBeamParticle) +
                    "].Description",
                "Incoming beam particle");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kTargetParticle) + "].Name",
                "TargetParticle");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kTargetParticle) +
                    "].Description",
                "The target particle in the hard scatter");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kPauliBlocked) + "].Name",
                "PauliBlocked");
  add_attribute(
      run_info,
      "NuHepMC.ParticleStatusInfo[" +
          std::to_string(ParticleStatus::kPauliBlocked) + "].Description",
      "Marks a particle as being Pauli blocked after the hard-scatter was "
      "simulated, generally this event should be skipped but included in "
      "total cross-section calculations.");

  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kStruckNucleon) + "].Name",
                "StruckNucleon");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kStruckNucleon) +
                    "].Description",
                "The nucleon involved in the hard scatter");

  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kUnderwentFSI) + "].Name",
                "UnderwentFSI");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kUnderwentFSI) +
                    "].Description",
                "This particle subsequently underwent FSI and should not be "
                "considered final state");

  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kSecondaryInteraction) +
                    "].Name",
                "SecondaryInteraction");
  add_attribute(run_info,
                "NuHepMC.ParticleStatusInfo[" +
                    std::to_string(ParticleStatus::kSecondaryInteraction) +
                    "].Description",
                "This particle subsequently underwent a secondary interaction "
                "in the detector");

  add_attribute(run_info, "NuHepMC.ParticleStatusIDs",
                std::vector<int>{ParticleStatus::kUndecayedPhysicalParticle,
                                 ParticleStatus::kDecayedParticle,
                                 ParticleStatus::kDocumentationLine,
                                 ParticleStatus::kIncomingBeamParticle,
                                 ParticleStatus::kTargetParticle,
                                 ParticleStatus::kPauliBlocked,
                                 ParticleStatus::kStruckNucleon,
                                 ParticleStatus::kUnderwentFSI,
                                 ParticleStatus::kSecondaryInteraction});

  // G.R.7 Event Weights
  run_info->set_weight_names({
      "CV",
  });

  // G.C.1 Signalling Followed Conventions
  add_attribute(run_info, "NuHepMC.Conventions",
                std::vector<std::string>{
                    "G.C.1",
                    "G.C.2",
                    "G.C.4"
                    "E.C.4",
                    "E.C.2",
                });

  // G.C.2 File Exposure (Standalone)
  add_attribute(run_info, "NuHepMC.Exposure.NEvents", nevents);

  // G.C.4 Flux-averaged Total Cross Section
  add_attribute(run_info, "NuHepMC.FluxAveragedTotalCrossSection",
                flux_averaged_total_cross_section);

  return run_info;
}

HepMC3::GenEvent ToGenEvent(NeutVect *nv, std::shared_ptr<HepMC3::GenRunInfo> gri) {

  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::CM);
  evt.set_run_info(gri);

  HepMC3::GenVertexPtr IAVertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  IAVertex->set_status(VertexStatus::kNuclearTargetVertex);

  int nuclear_PDG = nv->Ibound ? 1000000000 + nv->TargetZ * 10000 +
                                     (nv->TargetA + nv->TargetZ) * 10
                               : 1000010010;

  HepMC3::GenParticlePtr target_nucleus = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector{0, 0, 0, 0}, nuclear_PDG,
      ParticleStatus::kTargetParticle);

  IAVertex->add_particle_in(target_nucleus);

  // E.R.5
  HepMC3::GenVertexPtr primvertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  primvertex->set_status(VertexStatus::kPrimaryVertex);

  HepMC3::GenVertexPtr fsivertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  fsivertex->set_status(VertexStatus::kFSIVertex);

  NeutPart *pinfo = nullptr;
  int npart = nv->Npart();
  int nprimary = nv->Nprimary();
  for (int i = 0; i < npart; ++i) {

    pinfo = nv->PartInfo(i);

    bool isprim = i < nprimary;

    int NuHepPartStatus = 0;

    switch (pinfo->fStatus) {
    case -1: {
      NuHepPartStatus = (i == 0) ? ParticleStatus::kIncomingBeamParticle
                                 : ParticleStatus::kStruckNucleon;
      break;
    }
    case 0: {
      if (pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUndecayedPhysicalParticle;
      } else if ((std::abs(pinfo->fPID) == 12) ||
                 (std::abs(pinfo->fPID) == 14) ||
                 (std::abs(pinfo->fPID) == 16)) { // NC FS Neutrino
        NuHepPartStatus = ParticleStatus::kUndecayedPhysicalParticle;
      }
      break;
    }
    case 1: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kDecayedParticle;
      }
      break;
    }
    case 2: {
      if (pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUndecayedPhysicalParticle;
      } else if ((std::abs(pinfo->fPID) == 12) ||
                 (std::abs(pinfo->fPID) == 14) ||
                 (std::abs(pinfo->fPID) == 16)) { // NC FS Neutrino
        NuHepPartStatus = ParticleStatus::kUndecayedPhysicalParticle;
      }
      break;
    }
    case 3: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUnderwentFSI;
      }
      break;
    }
    case 4: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUnderwentFSI;
      }
      break;
    }
    case 5: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kPauliBlocked;
      }
      break;
    }
    case 6: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kSecondaryInteraction;
      }
      break;
    }
    case 7: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUnderwentFSI;
      }
      break;
    }
    case 8: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUnderwentFSI;
      }
      break;
    }
    case 9: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUnderwentFSI;
      }
      break;
    }
    case -3: { // absorbed pion
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = ParticleStatus::kUnderwentFSI;
      }
    }
    }

    if (!NuHepPartStatus) {
      std::cout << "[ERROR]: Failed to convert particle status for particle: "
                << i << std::endl;
      abort();
    }

    TLorentzVector fmom;
    fmom.SetXYZM(pinfo->fP.X(), pinfo->fP.Y(), pinfo->fP.Z(), pinfo->fMass);

    HepMC3::GenParticlePtr part = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{fmom.X(), fmom.Y(), fmom.Z(), fmom.E()}, pinfo->fPID,
        NuHepPartStatus);
    part->set_generated_mass(fmom.M());

    if (NuHepPartStatus == ParticleStatus::kIncomingBeamParticle) {
      primvertex->add_particle_in(part);
    } else if (NuHepPartStatus == ParticleStatus::kStruckNucleon) {
      primvertex->add_particle_in(part);
      IAVertex->add_particle_out(part);
    } else if (NuHepPartStatus == ParticleStatus::kUndecayedPhysicalParticle) {
      if (isprim) {
        primvertex->add_particle_out(part);
      } else {
        fsivertex->add_particle_out(part);
      }
    } else if ((NuHepPartStatus == ParticleStatus::kDecayedParticle) ||
               (NuHepPartStatus == ParticleStatus::kPauliBlocked) ||
               (NuHepPartStatus == ParticleStatus::kSecondaryInteraction)) {
      primvertex->add_particle_out(part);
    } else if (NuHepPartStatus == ParticleStatus::kUnderwentFSI) {
      primvertex->add_particle_out(part);
      fsivertex->add_particle_in(part);
    } else {
      std::cout << "[ERROR]: Failed to find vertex for particle: " << (i + 1)
                << std::endl;
      abort();
    }
  }

  if (!fsivertex->particles_in().size() && fsivertex->particles_out().size()) {
    // add a dummy documentation line connecting this vertex to the primary one
    HepMC3::GenParticlePtr doc_line = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0, 0, 0, 0}, 0, ParticleStatus::kDocumentationLine);
    primvertex->add_particle_out(doc_line);
    fsivertex->add_particle_in(doc_line);
  }

  evt.add_vertex(IAVertex);
  evt.add_vertex(primvertex);

  if (fsivertex->particles_in().size() || fsivertex->particles_out().size()) {
    evt.add_vertex(fsivertex);
  }

  // E.C.1
  evt.weight("CV") = 1;

  // E.C.4
  static double const cm2_to_pb = 1E36;

  // E.C.2
  add_attribute(evt, "TotXS", nv->Totcrs * 1E-38 * cm2_to_pb);

  // A semi random selection of properties that should exist in the last few
  // versions of NEUT
  add_attribute(evt, "TargetA", nv->TargetA);
  add_attribute(evt, "TargetZ", nv->TargetZ);
  add_attribute(evt, "TargetH", nv->TargetH);
  add_attribute(evt, "Ibound", nv->Ibound);
  add_attribute(evt, "VNuclIni", nv->VNuclIni);
  add_attribute(evt, "VNuclFin", nv->VNuclFin);
  add_attribute(evt, "PFSurf", nv->PFSurf);
  add_attribute(evt, "PFMax", nv->PFMax);
  add_attribute(evt, "QEModel", nv->QEModel);
  add_attribute(evt, "QEVForm", nv->QEVForm);
  add_attribute(evt, "RADcorr", nv->RADcorr);
  add_attribute(evt, "SPIModel", nv->SPIModel);
  add_attribute(evt, "COHModel", nv->COHModel);
  add_attribute(evt, "DISModel", nv->DISModel);

  // E.R.4
  add_attribute(evt, "LabPos", std::vector<double>{0, 0, 0, 0});

  add_attribute(evt, "ProcId", nv->Mode);

  return evt;
}