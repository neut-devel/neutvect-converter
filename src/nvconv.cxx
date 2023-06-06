#include "nvconv.h"

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

using namespace nvconverter;

// Lazy way of choosing the right attribute type via TMP
template <typename T> struct attr_traits {};

template <> struct attr_traits<int> {
  typedef HepMC3::IntAttribute type;
};

template <> struct attr_traits<double> {
  typedef HepMC3::DoubleAttribute type;
};

template <> struct attr_traits<std::string> {
  typedef HepMC3::StringAttribute type;
};

template <typename T>
void add_attribute(HepMC3::GenEvent &ge, std::string const &name,
                   T const &val) {
  ge.add_attribute(name, std::make_shared<typename attr_traits<T>::type>(val));
}

static const double GeV = 1E-3;
static const double MeV = 1;

HepMC3::GenEvent ToGenEvent(NeutVect *nv,
                            double flux_averaged_total_cross_section) {

  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::CM);

  HepMC3::GenVertexPtr hard_scatter_vtx =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{0, 0, 0, 0});

  HepMC3::GenVertexPtr fsi_vtx =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector{0, 0, 0, 0});

  NeutPart *pinfo = nullptr;
  int npart = nv->Npart();
  int nprimary = nv->Nprimary();
  for (int i = 0; i < npart; ++i) {
    pinfo = nv->PartInfo(i);

    HepMC3::GenParticlePtr part = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector{pinfo->fP.X() * GeV, pinfo->fP.Y() * GeV,
                           pinfo->fP.Z() * GeV, pinfo->fP.E() * GeV},
        pinfo->fPID, 0);

    part->set_generated_mass(pinfo->fMass * GeV);

    if (pinfo->fStatus == -1) { // Inital State

      part->set_status(pinfo->fPID > 10000 ? kTgtNucleusStatus : kISStatus);
      hard_scatter_vtx->add_particle_in(part);

    } else if (pinfo->fStatus == 0) {

      if (pinfo->fIsAlive) {
        part->set_status(kFSStatus);
      } else {
        part->set_status(kOther);
      }

      if (i < nprimary) {
        hard_scatter_vtx->add_particle_out(part);
      }
      fsi_vtx->add_particle_out(part);

    } else {

      switch (pinfo->fStatus) {
      case 1: {
        part->set_status(kFSIDecayed);
        break;
      }
      case 3: {
        part->set_status(kFSIAbsorpbed);
        break;
      }
      case 4: {
        part->set_status(kFSIChargeExchanged);
        break;
      }
      case 5: {
        part->set_status(kPauliBlocked);
        break;
      }
      }

      if (i < nprimary) {
        hard_scatter_vtx->add_particle_out(part);
      }
      fsi_vtx->add_particle_out(part);
    }
  }

  evt.add_vertex(hard_scatter_vtx);
  evt.add_vertex(fsi_vtx);

  add_attribute(evt, "ProcId", nv->Mode);

  return evt;
}