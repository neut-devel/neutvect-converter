#include "nvconv.h"

#include "NuHepMC/WriterUtils.hxx"

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#ifdef NEUTCONV_DEBUG
#include "HepMC3/Print.h"
#endif

#include <memory>
#include <set>
#include <utility>

namespace nvconv {

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

void AddNEUTPassthrough(HepMC3::GenEvent &evt,
                        std::vector<HepMC3::GenParticlePtr> &parts,
                        NeutVect *nv) {

  NuHepMC::add_attribute(evt, "NEUT.TargetA", nv->TargetA);
  NuHepMC::add_attribute(evt, "NEUT.TargetZ", nv->TargetZ);
  NuHepMC::add_attribute(evt, "NEUT.TargetH", nv->TargetH);
  NuHepMC::add_attribute(evt, "NEUT.Ibound", nv->Ibound);
  NuHepMC::add_attribute(evt, "NEUT.VNuclIni", nv->VNuclIni);
  NuHepMC::add_attribute(evt, "NEUT.VNuclFin", nv->VNuclFin);
  NuHepMC::add_attribute(evt, "NEUT.PFSurf", nv->PFSurf);
  NuHepMC::add_attribute(evt, "NEUT.PFMax", nv->PFMax);
  NuHepMC::add_attribute(evt, "NEUT.QEModel", nv->QEModel);
  NuHepMC::add_attribute(evt, "NEUT.QEVForm", nv->QEVForm);
  NuHepMC::add_attribute(evt, "NEUT.RADcorr", nv->RADcorr);
  NuHepMC::add_attribute(evt, "NEUT.SPIModel", nv->SPIModel);
  NuHepMC::add_attribute(evt, "NEUT.COHModel", nv->COHModel);
  NuHepMC::add_attribute(evt, "NEUT.DISModel", nv->DISModel);
  NuHepMC::add_attribute(evt, "NEUT.Mode", nv->Mode);

  NeutPart *pinfo = nullptr;
  int npart = nv->Npart();
  int nprimary = nv->Nprimary();

  NuHepMC::add_attribute(evt, "NEUT.npart", npart);
  NuHepMC::add_attribute(evt, "NEUT.nprimary", nprimary);

  for (int i = 0; i < npart; ++i) {
    pinfo = nv->PartInfo(i);
    if (parts[i]->in_event()) {
      NuHepMC::add_attribute(parts[i], "NEUT.i", i);
      NuHepMC::add_attribute(parts[i], "NEUT.fStatus", pinfo->fStatus);
      NuHepMC::add_attribute(parts[i], "NEUT.fIsAlive", pinfo->fIsAlive);
    }
  }
}

std::shared_ptr<HepMC3::GenRunInfo>
BuildRunInfo(int nevents, double flux_averaged_total_cross_section,
             std::unique_ptr<TH1> &flux_hist, bool &isMonoE, int beam_pid,
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

  NuHepMC::GR9::WriteVertexStatusIDDefinitions(run_info, VertexStatuses);

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
  NuHepMC::GR10::WriteParticleStatusIDDefinitions(run_info, ParticleStatuses);

  // G.R.7 Event Weights
  NuHepMC::GR7::SetWeightNames(run_info, {
                                             "CV",
                                         });

  // G.R.4 Signalling Followed Conventions
  std::vector<std::string> conventions = {
      "G.C.2", "E.C.1", "E.C.2", "V.C.1", "P.C.1", "P.C.2",
  };

  // G.R.6 Cross Section Units and Target Scaling
  NuHepMC::GR6::SetCrossSectionUnits(run_info, "pb", "PerNucleon");

  // G.C.2 Flux-averaged Total Cross Section
  NuHepMC::GC2::SetFluxAveragedTotalXSec(run_info,
                                         flux_averaged_total_cross_section);

  if (flux_hist || isMonoE) {
    conventions.push_back("G.C.7");

    if (isMonoE) {

      NuHepMC::GC4::WriteBeamUnits(run_info, "MEV");
      NuHepMC::GC4::SetMonoEnergeticBeamType(run_info);

    } else {

      NuHepMC::GC4::SetHistogramBeamType(run_info);

      std::vector<double> bin_edges;
      std::vector<double> bin_content;

      bin_edges.push_back(flux_hist->GetXaxis()->GetBinLowEdge(1) *
                          flux_to_MeV);
      for (int i = 0; i < flux_hist->GetXaxis()->GetNbins(); ++i) {
        bin_edges.push_back(flux_hist->GetXaxis()->GetBinUpEdge(i + 1) *
                            flux_to_MeV);
        bin_content.push_back(flux_hist->GetBinContent(i + 1));
      }

      NuHepMC::GC4::WriteBeamEnergyHistogram(run_info, beam_pid, bin_edges,
                                             bin_content, false);
      NuHepMC::GC4::WriteBeamUnits(
          run_info, "MEV", std::string(flux_hist->GetYaxis()->GetTitle()));
    }
  }

  // G.R.4 Signalling Followed Conventions
  NuHepMC::GR4::SetConventions(run_info, conventions);

  return run_info;
}

std::shared_ptr<HepMC3::GenEvent>
ToGenEvent(NeutVect *nv, std::shared_ptr<HepMC3::GenRunInfo> gri) {

#ifdef NEUTCONV_DEBUG
  std::cout << ">>>>>>>>>>>>>>>ToGenEvent" << std::endl;
  nv->Dump();
#endif

  auto evt =
      std::make_shared<HepMC3::GenEvent>(HepMC3::Units::MEV, HepMC3::Units::CM);
  evt->set_run_info(gri);

  HepMC3::GenVertexPtr IAVertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  IAVertex->set_status(NuHepMC::VertexStatus::NucleonSeparation);

  bool isbound = nv->Ibound;

  std::set<int> not_bound_modes = {16, 15, -16, -15, 36, 35, -36, -35};
  if (!isbound && (not_bound_modes.count(nv->Mode))) {
    // Correct confusing 'isbound == false' for certain modes
    isbound = true;
  }

  int nuclear_PDG = 1000000000 + nv->TargetZ * 10000 + nv->TargetA * 10;
  int nuclear_remnant_PDG = nuclear_PDG;
  HepMC3::GenParticlePtr nuclear_remnant_internal = nullptr;
  HepMC3::GenParticlePtr nuclear_remnant_external = nullptr;
  if (isbound) {

    HepMC3::GenParticlePtr target_nucleus =
        std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector{0, 0, 0, 0},
                                              nuclear_PDG,
                                              NuHepMC::ParticleStatus::Target);

    nuclear_remnant_internal = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0, 0, 0, 0}, NuHepMC::ParticleNumber::NuclearRemnant,
        NuHepMC::ParticleStatus::DocumentationLine);

    nuclear_remnant_external = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{0, 0, 0, 0}, NuHepMC::ParticleNumber::NuclearRemnant,
        NuHepMC::ParticleStatus::UndecayedPhysical);

    IAVertex->add_particle_in(target_nucleus);
    IAVertex->add_particle_out(nuclear_remnant_internal);
  }

  // E.R.5
  HepMC3::GenVertexPtr primvertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  primvertex->set_status(NuHepMC::VertexStatus::Primary);

  HepMC3::GenVertexPtr fsivertex =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{});
  fsivertex->set_status(NuHepMC::VertexStatus::FSISummary);
  fsivertex->add_particle_in(nuclear_remnant_internal);
  fsivertex->add_particle_out(nuclear_remnant_external);

  NeutPart *pinfo = nullptr;
  int npart = nv->Npart();
  int nprimary = nv->Nprimary();

  // need to keep this stack so that we can add metadata attributes after we
  // have added them to the event.
  std::vector<HepMC3::GenParticlePtr> parts;

  for (int p_it = 0; p_it < npart; ++p_it) {
    pinfo = nv->PartInfo(p_it);

    bool isprim = p_it < nprimary;

    int NuHepPartStatus = 0;

    switch (pinfo->fStatus) {
    case -1: {
      NuHepPartStatus = (p_it == 0) ? NuHepMC::ParticleStatus::IncomingBeam
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

    if ((std::abs(nv->Mode) == 15) ||
        (std::abs(nv->Mode) == 35)) { // special case for diffractive
      switch (p_it) {
      case 0: {
        NuHepPartStatus = NuHepMC::ParticleStatus::IncomingBeam;
        break;
      }
      case 1: {
        NuHepPartStatus = NuHepMC::ParticleStatus::StruckNucleon;
        break;
      }
      default: {
        NuHepPartStatus = NuHepMC::ParticleStatus::UndecayedPhysical;
        break;
      }
      }
    }

    if (!NuHepPartStatus) {

      std::stringstream ss;
      ss << "[ERROR]: Failed to convert particle status for particle: " << p_it
         << "\n";

      for (int p_it = 0; p_it < npart; ++p_it) {
        auto pinfo = nv->PartInfo(p_it);
        ss << "p[" << p_it << "]- pid: " << pinfo->fPID
           << ", prim: " << (p_it < nprimary) << ", status: " << pinfo->fStatus
           << ", alive: " << pinfo->fIsAlive << "\n";
      }

      std::cout << ss.str() << std::endl;
      nv->Dump();
      throw ss.str();
    }

    TLorentzVector fmom;
    fmom.SetXYZM(pinfo->fP.X(), pinfo->fP.Y(), pinfo->fP.Z(), pinfo->fMass);

    HepMC3::GenParticlePtr part = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{fmom.X(), fmom.Y(), fmom.Z(), fmom.E()}, pinfo->fPID,
        NuHepPartStatus);
    parts.push_back(part);
    part->set_generated_mass(fmom.M());

#ifdef NEUTCONV_DEBUG
    std::cout << "Processing NEUT particle(" << p_it << "/" << npart
              << ", nprim:" << nprimary << ") pid: " << pinfo->fPID
              << ", status: " << pinfo->fStatus
              << " isalive: " << pinfo->fIsAlive << std::endl;
    std::cout << "\t->NuHepPartStatus: " << NuHepPartStatus << std::endl;
#endif

    if (NuHepPartStatus == NuHepMC::ParticleStatus::IncomingBeam) {
      primvertex->add_particle_in(part);
#ifdef NEUTCONV_DEBUG
      std::cout << "\t->Added as /in/ to primvertex" << std::endl;
#endif
    } else if (NuHepPartStatus == NuHepMC::ParticleStatus::StruckNucleon) {
      if (isbound) {
        IAVertex->add_particle_out(part);

        if (part->pid() == 2212) {
          nuclear_remnant_PDG -= (1 * 10000 + 1 * 10);
        } else {
          nuclear_remnant_PDG -= (0 * 10000 + 1 * 10);
        }
      } else { // use stuck nucleon as target for unbound interactions
        part->set_status(NuHepMC::ParticleStatus::Target);
      }
      primvertex->add_particle_in(part);

#ifdef NEUTCONV_DEBUG
      if (isbound) {
        std::cout << "\t->Added as /out/ from IAvertex" << std::endl;
      }
      std::cout << "\t->Added as /in/ to primvertex" << std::endl;
#endif
    } else if (NuHepPartStatus == NuHepMC::ParticleStatus::UndecayedPhysical) {
      if (isprim) {
        auto part_copy = std::make_shared<HepMC3::GenParticle>(part->data());
        part_copy->set_status(NuHepMC::ParticleStatus::DocumentationLine);
        primvertex->add_particle_out(part_copy);
        if (isbound) {
          fsivertex->add_particle_in(part_copy);
        }
#ifdef NEUTCONV_DEBUG
        std::cout << "\t->Copied particle with status: "
                  << NuHepMC::ParticleStatus::DocumentationLine << std::endl;
        std::cout << "\t\t->Added as /out/ from primvertex" << std::endl;
        if (isbound) {
          std::cout << "\t\t->Added as /in/ to fsivertex" << std::endl;
        }
#endif
      }
      if (isbound) {
        fsivertex->add_particle_out(part);
#ifdef NEUTCONV_DEBUG
        std::cout << "\t->Added as /out/ from fsivertex" << std::endl;
#endif
      }
    } else if ((NuHepPartStatus == NuHepMC::ParticleStatus::DecayedPhysical) ||
               (NuHepPartStatus ==
                NuHepMC::ParticleStatus::NEUT::PauliBlocked) ||
               (NuHepPartStatus ==
                NuHepMC::ParticleStatus::NEUT::SecondaryInteraction)) {
      primvertex->add_particle_out(part);
      fsivertex->add_particle_in(part);
#ifdef NEUTCONV_DEBUG
      std::cout << "\t->Added as /out/ from primvertex" << std::endl;
      std::cout << "\t->Added as /in/ to fsivertex" << std::endl;
#endif
    } else if (NuHepPartStatus == NuHepMC::ParticleStatus::NEUT::UnderwentFSI) {
      primvertex->add_particle_out(part);
      fsivertex->add_particle_in(part);
#ifdef NEUTCONV_DEBUG
      std::cout << "\t->Added as /out/ from primvertex" << std::endl;
      std::cout << "\t->Added as /in/ to fsivertex" << std::endl;
#endif
    } else {
      std::stringstream ss;
      ss << "[ERROR]: Failed to find vertex for particle: " << (p_it + 1);
      std::cout << ss.str() << std::endl;
      throw ss.str();
    }
  }

  if (isbound && !fsivertex->particles_in().size()) {
    throw std::runtime_error(
        "neutvect-converter: [ERROR]: fsivertex had no incoming particles.");
  }

  if (isbound) {
    evt->add_vertex(IAVertex);
  }
  evt->add_vertex(primvertex);
  if (isbound) {
    evt->add_vertex(fsivertex);
  }

  // E.C.1
  evt->weight("CV") = 1;

  // E.C.4
  static double const cm2_to_pb = 1E36;

  // E.C.2
  NuHepMC::EC2::SetTotalCrossSection(*evt, nv->Totcrs * 1E-38 * cm2_to_pb);
  // E.R.5
  NuHepMC::ER5::SetLabPosition(*evt, std::vector<double>{0, 0, 0, 0});

  NuHepMC::ER3::SetProcessID(*evt, GetEC1Channel(nv->Mode));

  if (isbound) {
    NuHepMC::PC2::SetRemnantNucleusParticleNumber(
        nuclear_remnant_internal, (nuclear_remnant_PDG / 10000) % 1000,
        (nuclear_remnant_PDG / 10) % 1000);
  }

  AddNEUTPassthrough(*evt, parts, nv);

#ifdef NEUTCONV_DEBUG
  HepMC3::Print::listing(*evt);
  std::cout << "<<<<<<<<<<<<<<<ToGenEvent" << std::endl;
#endif

  return evt;
}
} // namespace nvconv