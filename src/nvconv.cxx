#include "nvconv.h"

#include "NuHepMC/WriterUtils.hxx"

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include <set>
#include <utility>

static const double GeV = 1E-3;
static const double MeV = 1;

const std::map<int, std::pair<std::string, int>> ChannelNameIndexModeMapping{
    {1, {"CC_QE_nu", 200}},

    {11, {"CC_RES_ppi+_nu", 400}},
    {12, {"CC_RES_ppi0_nu", 401}},
    {13, {"CC_RES_npi+_nu", 402}},

    {21, {"CC_multi_pi_nu", 500}},

    {31, {"NC_RES_npi0_nu", 450}},
    {32, {"NC_RES_ppi0_nu", 451}},
    {33, {"NC_RES_ppi-_nu", 452}},
    {34, {"NC_RES_npi+_nu", 452}},

    {41, {"NC_multi_pi_nu", 550}},

    {51, {"NC_elastic_p_nu", 250}},
    {52, {"NC_elastic_n_nu", 251}},

    {16, {"CC_COH_nu", 100}},

    {36, {"NC_COH_nu", 150}},

    {22, {"CC_eta_nu", 410}},

    {42, {"NC_eta_n_nu", 460}},
    {43, {"NC_eta_p_nu", 461}},

    {23, {"CC_kaon_nu", 411}},

    {44, {"NC_kaon_n_nu", 462}},
    {45, {"NC_kaon_p_nu", 463}},

    {26, {"CC_DIS_nu", 600}},

    {46, {"NC_DIS_nu", 601}},

    {17, {"CC_1gamma_nu", 412}},

    {38, {"NC_1gamma_n_nu", 464}},
    {39, {"NC_1gamma_p_nu", 465}},

    {2, {"CC_2p2h_nu", 300}},

    {15, {"CC_DIF_nu", 110}},

    {35, {"NC_DIF_nu", 160}},

    {-1, {"CC_QE_proton_nubar", 225}},

    {-11, {"CC_RES_npi-_nubar", 425}},
    {-12, {"CC_RES_ppi0_nubar", 426}},
    {-13, {"CC_RES_ppi-_nubar", 427}},

    {-21, {"CC_multi_pi_nubar", 525}},

    {-31, {"NC_RES_npi0_nubar", 475}},
    {-32, {"NC_RES_ppi0_nubar", 476}},
    {-33, {"NC_RES_ppi-_nubar", 478}},
    {-34, {"NC_RES_npi+_nubar", 479}},

    {-41, {"NC_multi_pi_nubar", 575}},

    {-51, {"NC_elastic_p_nubar", 275}},
    {-52, {"NC_elastic_n_nubar", 276}},

    {-16, {"CC_COH_nubar", 125}},

    {-36, {"NC_COH_nubar", 175}},

    {-22, {"CC_eta_nubar", 435}},

    {-42, {"NC_eta_n_nubar", 485}},
    {-43, {"NC_eta_p_nubar", 486}},

    {-23, {"CC_kaon_nubar", 436}},

    {-44, {"NC_kaon_n_nubar", 487}},
    {-45, {"NC_kaon_p_nubar", 488}},

    {-26, {"CC_DIS_nubar", 625}},

    {-46, {"NC_DIS_nubar", 675}},

    {-17, {"CC_1gamma_nubar", 437}},

    {-38, {"NC_1gamma_n_nubar", 489}},
    {-39, {"NC_1gamma_p_nubar", 490}},

    {-2, {"CC_2p2h_nubar", 325}},

    {-15, {"CC_DIF_nubar", 135}},

    {-35, {"NC_DIF_nubar", 185}},

};

int GetEC1Channel(int neutmode) {
  if (!ChannelNameIndexModeMapping.count(neutmode)) {
    std::cout << "[ERROR]: neutmode: " << neutmode << " unaccounted for."
              << std::endl;
    throw neutmode;
  }
  return ChannelNameIndexModeMapping.at(neutmode).second;
}

std::shared_ptr<HepMC3::GenRunInfo>
BuildRunInfo(int nevents, double flux_averaged_total_cross_section,
             TH1 *&flux_histo, bool &isMonoE, int beam_pid,
             double flux_to_MeV) {

  // G.R.1 Valid GenRunInfo
  auto run_info = std::make_shared<HepMC3::GenRunInfo>();

  // G.R.2 NuHepMC Version
  NuHepMC::GR2::WriteVersion(run_info);

  // G.R.3 Generator Identification
  run_info->tools().emplace_back(HepMC3::GenRunInfo::ToolInfo{
      "NEUT", NEUT_VERSION_STR,
      "https://doi.org/10.1140/epjs/s11734-021-00287-7"});
  run_info->tools().emplace_back(
      HepMC3::GenRunInfo::ToolInfo{"neutvect-converter", PROJECT_VERSION_STR,
                                   "github.com/neut-dev/neutvect-converter"});

  // G.R.4 Process Metadata

  std::vector<int> pids;
  for (auto const &neutchan : ChannelNameIndexModeMapping) {
    pids.push_back(neutchan.second.second);
    NuHepMC::add_attribute(run_info,
                           "NuHepMC.ProcessInfo[" +
                               std::to_string(neutchan.second.second) +
                               "].Name",
                           neutchan.second.first);
    NuHepMC::add_attribute(
        run_info,
        "NuHepMC.ProcessInfo[" + std::to_string(neutchan.second.second) +
            "].Description",
        std::string("neutmode=") + std::to_string(neutchan.first));
  }

  NuHepMC::add_attribute(run_info, "NuHepMC.ProcessIDs", pids);

  // G.R.5 Vertex Status Metadata
  NuHepMC::StatusCodeDescriptors VertexStatuses = {
      {NuHepMC::VertexStatus::Primary,
       {"PrimaryVertex", "The neutrino hard-scatter vertex"}},
      {NuHepMC::VertexStatus::FSISummary,
       {"FSIVertex", "A single vertex representing the cascade"}},
      {NuHepMC::VertexStatus::NucleonSeparation,
       {"NucleonSeparationVertex",
        "Impulse approximation vertex that represents the separation of the "
        "single target nucleon from the target nucleus ground state."}},
  };

  NuHepMC::GR5::WriteVertexStatusIDDefinitions(run_info, VertexStatuses);

  NuHepMC::StatusCodeDescriptors ParticleStatuses = {
      {NuHepMC::ParticleStatus::UndecayedPhysical,
       {"UndecayedPhysical",
        "Physical final state particles produced by this simulation"}},
      {NuHepMC::ParticleStatus::DecayedPhysical,
       {"DecayedPhysical", "Particle was decayed by the simulation"}},
      {NuHepMC::ParticleStatus::DocumentationLine,
       {"DocumentationLine",
        "Documentation line, not considered a real particle"}},
      {NuHepMC::ParticleStatus::IncomingBeam,
       {"IncomingBeam", "Incoming beam particle"}},
      {NuHepMC::ParticleStatus::Target,
       {"TargetParticle", "The target particle in the hard scatter"}},
      {NuHepMC::ParticleStatus::StruckNucleon,
       {"StruckNucleon", "The nucleon involved in the hard scatter"}},
      {NuHepMC::ParticleStatus::NEUT::PauliBlocked,
       {"PauliBlocked",
        "Marks a particle as being Pauli blocked after the hard-scatter was "
        "simulated, generally this event should be skipped but included in "
        "total "
        "cross-section calculations."}},
      {NuHepMC::ParticleStatus::NEUT::UnderwentFSI,
       {"UnderwentFSI", "This particle subsequently underwent FSI and should "
                        "not be considered part of the final state"}},
      {NuHepMC::ParticleStatus::NEUT::SecondaryInteraction,
       {"SecondaryInteraction", "This particle subsequently underwent a "
                                "secondary interaction in the detector"}}};
  NuHepMC::GR6::WriteParticleStatusIDDefinitions(run_info, ParticleStatuses);

  // G.R.7 Event Weights
  NuHepMC::GR7::SetWeightNames(run_info, {
                                             "CV",
                                         });

  // G.C.1 Signalling Followed Conventions
  std::vector<std::string> conventions = {
      "G.C.1", "G.C.2", "G.C.4", "G.C.5", "E.C.1", "E.C.2", "E.C.3",
  };

  // G.C.2 File Exposure (Standalone)
  NuHepMC::GC2::SetExposureNEvents(run_info, nevents);

  // G.C.4 Cross Section Units and Target Scaling
  NuHepMC::GC4::SetCrossSectionUnits(run_info, "pb", "PerTargetMolecule");

  // G.C.5 Flux-averaged Total Cross Section
  NuHepMC::GC5::SetFluxAveragedTotalXSec(run_info,
                                         flux_averaged_total_cross_section);

  if (flux_histo || isMonoE) {
    conventions.push_back("G.C.7");

    if (isMonoE) {
      NuHepMC::add_attribute(
          run_info, "NuHepMC.Beam[" + std::to_string(beam_pid) + "].Type",
          "MonoEnergetic");

      NuHepMC::add_attribute(
          run_info, "NuHepMC.Beam[" + std::to_string(beam_pid) + "].EnergyUnit",
          "MEV");

    } else {
      NuHepMC::add_attribute(
          run_info, "NuHepMC.Beam[" + std::to_string(beam_pid) + "].Type",
          "Histogram");

      std::vector<double> bin_edges;
      std::vector<double> bin_content;

      bin_edges.push_back(flux_histo->GetXaxis()->GetBinLowEdge(1) *
                          flux_to_MeV);
      for (int i = 0; i < flux_histo->GetXaxis()->GetNbins(); ++i) {
        bin_edges.push_back(flux_histo->GetXaxis()->GetBinUpEdge(i + 1) *
                            flux_to_MeV);
        bin_content.push_back(flux_histo->GetBinContent(i + 1));
      }

      NuHepMC::add_attribute(run_info,
                             "NuHepMC.Beam[" + std::to_string(beam_pid) +
                                 "].Histogram.BinEdges",
                             bin_edges);

      NuHepMC::add_attribute(
          run_info, "NuHepMC.Beam[" + std::to_string(beam_pid) + "].EnergyUnit",
          "MEV");

      if (std::string(flux_histo->GetYaxis()->GetTitle()).size()) {
        NuHepMC::add_attribute(
            run_info, "NuHepMC.Beam[" + std::to_string(beam_pid) + "].RateUnit",
            flux_histo->GetYaxis()->GetTitle());
      }

      NuHepMC::add_attribute<bool>(run_info,
                                   "NuHepMC.Beam[" + std::to_string(beam_pid) +
                                       "].Histogram.ContentIsPerWidth",
                                   false);

      NuHepMC::add_attribute(run_info,
                             "NuHepMC.Beam[" + std::to_string(beam_pid) +
                                 "].Histogram.BinContent",
                             bin_content);
    }
  }

  // G.C.1 Signalling Followed Conventions
  NuHepMC::GC1::SetConventions(run_info, conventions);

  return run_info;
}

HepMC3::GenEvent ToGenEvent(NeutVect *nv,
                            std::shared_ptr<HepMC3::GenRunInfo> gri) {

  HepMC3::GenEvent evt(HepMC3::Units::MEV, HepMC3::Units::CM);
  evt.set_run_info(gri);

  HepMC3::GenVertexPtr IAVertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  IAVertex->set_status(NuHepMC::VertexStatus::NucleonSeparation);

  bool isbound = nv->Ibound;

  std::set<int> not_bound_modes = {16, -16, 36, -36};
  if (!isbound && (not_bound_modes.count(nv->Mode))) {
    // Correct confusing 'isbound == false' for certain modes
    isbound = true;
  }

  int nuclear_PDG = isbound
                        ? 1000000000 + nv->TargetZ * 10000 + nv->TargetA * 10
                        : 1000010010;

  HepMC3::GenParticlePtr target_nucleus = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector{0, 0, 0, 0}, nuclear_PDG,
      NuHepMC::ParticleStatus::Target);

  IAVertex->add_particle_in(target_nucleus);

  // E.R.5
  HepMC3::GenVertexPtr primvertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  primvertex->set_status(NuHepMC::VertexStatus::Primary);

  HepMC3::GenVertexPtr fsivertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  fsivertex->set_status(NuHepMC::VertexStatus::FSISummary);

  NeutPart *pinfo = nullptr;
  int npart = nv->Npart();
  int nprimary = nv->Nprimary();
  for (int i = 0; i < npart; ++i) {
    pinfo = nv->PartInfo(i);

    bool isprim = i < nprimary;

    int NuHepPartStatus = 0;

    switch (pinfo->fStatus) {
    case -1: {
      NuHepPartStatus = (i == 0) ? NuHepMC::ParticleStatus::IncomingBeam
                                 : NuHepMC::ParticleStatus::StruckNucleon;
      break;
    }
    case 0: {
      if (pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::UndecayedPhysical;
      } else if ((std::abs(pinfo->fPID) == 12) ||
                 (std::abs(pinfo->fPID) == 14) ||
                 (std::abs(pinfo->fPID) == 16)) { // NC FS Neutrino
        NuHepPartStatus = NuHepMC::ParticleStatus::UndecayedPhysical;
      }
      break;
    }
    case 1: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::DecayedPhysical;
      }
      break;
    }
    case 2: {
      if (pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::UndecayedPhysical;
      } else if ((std::abs(pinfo->fPID) == 12) ||
                 (std::abs(pinfo->fPID) == 14) ||
                 (std::abs(pinfo->fPID) == 16)) { // NC FS Neutrino
        NuHepPartStatus = NuHepMC::ParticleStatus::UndecayedPhysical;
      }
      break;
    }
    case 3: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::UnderwentFSI;
      }
      break;
    }
    case 4: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::UnderwentFSI;
      }
      break;
    }
    case 5: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::PauliBlocked;
      }
      break;
    }
    case 6: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::SecondaryInteraction;
      }
      break;
    }
    case 7: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::UnderwentFSI;
      }
      break;
    }
    case 8: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::UnderwentFSI;
      }
      break;
    }
    case 9: {
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::UnderwentFSI;
      }
      break;
    }
    case -3: { // absorbed pion
      if (!pinfo->fIsAlive) {
        NuHepPartStatus = NuHepMC::ParticleStatus::NEUT::UnderwentFSI;
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

    if (NuHepPartStatus == NuHepMC::ParticleStatus::IncomingBeam) {
      primvertex->add_particle_in(part);
    } else if (NuHepPartStatus == NuHepMC::ParticleStatus::StruckNucleon) {
      primvertex->add_particle_in(part);
      IAVertex->add_particle_out(part);
    } else if (NuHepPartStatus == NuHepMC::ParticleStatus::UndecayedPhysical) {
      if (isprim) {
        primvertex->add_particle_out(part);
      } else {
        fsivertex->add_particle_out(part);
      }
    } else if ((NuHepPartStatus == NuHepMC::ParticleStatus::DecayedPhysical) ||
               (NuHepPartStatus ==
                NuHepMC::ParticleStatus::NEUT::PauliBlocked) ||
               (NuHepPartStatus ==
                NuHepMC::ParticleStatus::NEUT::SecondaryInteraction)) {
      primvertex->add_particle_out(part);
    } else if (NuHepPartStatus == NuHepMC::ParticleStatus::NEUT::UnderwentFSI) {
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
        HepMC3::FourVector{0, 0, 0, 0}, 0,
        NuHepMC::ParticleStatus::DocumentationLine);
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
  NuHepMC::EC2::SetTotalCrossSection(evt, nv->Totcrs * 1E-38 * cm2_to_pb *
                                              (nv->TargetA + nv->TargetH));

  // A semi random selection of properties that should exist in the last few
  // versions of NEUT
  NuHepMC::add_attribute(evt, "TargetA", nv->TargetA);
  NuHepMC::add_attribute(evt, "TargetZ", nv->TargetZ);
  NuHepMC::add_attribute(evt, "TargetH", nv->TargetH);
  NuHepMC::add_attribute(evt, "Ibound", nv->Ibound);
  NuHepMC::add_attribute(evt, "VNuclIni", nv->VNuclIni);
  NuHepMC::add_attribute(evt, "VNuclFin", nv->VNuclFin);
  NuHepMC::add_attribute(evt, "PFSurf", nv->PFSurf);
  NuHepMC::add_attribute(evt, "PFMax", nv->PFMax);
  NuHepMC::add_attribute(evt, "QEModel", nv->QEModel);
  NuHepMC::add_attribute(evt, "QEVForm", nv->QEVForm);
  NuHepMC::add_attribute(evt, "RADcorr", nv->RADcorr);
  NuHepMC::add_attribute(evt, "SPIModel", nv->SPIModel);
  NuHepMC::add_attribute(evt, "COHModel", nv->COHModel);
  NuHepMC::add_attribute(evt, "DISModel", nv->DISModel);

  // E.R.4
  NuHepMC::ER4::SetLabPosition(evt, std::vector<double>{0, 0, 0, 0});

  NuHepMC::ER2::SetProcessID(evt, GetEC1Channel(nv->Mode));

  return evt;
}