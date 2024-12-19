/**
 * @file utilities_pi02024.h
 * @brief Header file for definitions of utility functions for supporting
 * analysis variables and cuts.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author lkashur@colostate.edu
 */
#ifndef UTILITIES_PI02024_H
#define UTILITIES_PI02024_H

#include <vector>
#include <TVector3.h>
#include "include/cuts.h"

#define MIN_PHOTON_ENERGY 25.0
#define MIN_MUON_ENERGY 143.425
#define MIN_PION_ENERGY 25.0

struct truth_inter {
  int num_primary_muons;
  int num_primary_muons_thresh;
  int num_primary_pions;
  int num_primary_pions_thresh;
  int num_primary_pi0s;
  int num_primary_pi0s_thresh;
  int num_nonprimary_pi0s;
  double transverse_momentum_mag;
  bool is_fiducial;
  bool has_contained_tracks;
  bool is_neutrino;
  bool is_cc;
  double muon_momentum_mag;
  double muon_beam_costheta;
  double pi0_leading_photon_energy;
  double pi0_leading_photon_conv_dist;
  double pi0_subleading_photon_energy;
  double pi0_subleading_photon_conv_dist;
  double pi0_costheta;
  double pi0_mass;
  double pi0_momentum_mag;
  double pi0_beam_costheta;
};

struct reco_inter {
  double transverse_momentum_mag;
  double muon_momentum_mag;
  double muon_beam_costheta;
  double pi0_leading_photon_energy;
  double pi0_leading_photon_cosphi;
  double pi0_leading_photon_conv_dist;
  double pi0_subleading_photon_energy;
  size_t pi0_subleading_photon_cosphi;
  double pi0_subleading_photon_conv_dist;
  double pi0_costheta;
  double pi0_mass;
  double pi0_momentum_mag;
  double pi0_beam_costheta;
};



/**
 * @namespace utilities_pi02024
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities_pi02024
{
    /**
     * @brief Check if the particle meets final state signal requirements.
     * @details Particle must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), pions
     * an energy of at least 25 MeV, and photons an energy of at least 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
        bool final_state_signal(const T & p)
        {
	  bool passes(false);
	  if(p.is_primary)
	  {
	    double energy(p.pid > 1 ? p.csda_ke : p.calo_ke); // only valid for contained analysis  
	    if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
			     energy = p.ke;

	    if(p.pid == 0 && energy >= MIN_PHOTON_ENERGY) passes = true; // Photons
	    if(p.pid == 1) passes = true; // Electrons
	    if(p.pid == 2 && energy >= MIN_MUON_ENERGY) passes = true; // Muons
	    if(p.pid == 3 && energy >= MIN_PION_ENERGY) passes = true; // Pions
	    if(p.pid == 4) passes = true; // Protons                                        
	    if(p.pid == 5) passes = true; // Kaons
	  }
          return passes;
	}

    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the ineraction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
	  std::vector<uint32_t> counts(6, 0);
	  for(auto &p : obj.particles)
          {
            if(final_state_signal(p))
              ++counts[p.pid];
          }
	  return counts;
	}

    /**
     * @brief Structure for true interaction information.
     * @details This structure is needed because true neutral
     * pions are not kept track of in SPINE labeling.  Instead,
     * photons belonging to neutral pions are grouped in pairs.
     * @tparam T the type of interaction (true).
     * @param obj the interaction to select on.
     * @return a truth_inter structure.
     * @note This structure is intented to be used for the pi02024 analysis. 
     */
    template<class T> 
      truth_inter truth_interaction_info(const T & obj)
      {
	// Initialize struct
	truth_inter s;
	  
	// Initialize relevant TVector3s
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
	TVector3 beamdir(0, 0, 1); // BNB
	  
	// Initialize output variables
	int primary_muon_count(0);
	int primary_muon_count_thresh(0);
	int primary_pion_count(0);
	int primary_pion_count_thresh(0);
	int primary_pi0_count(0);
	int primary_pi0_count_thresh(0);
	int nonprimary_pi0_count(0);
	double pT0(0), pT1(0), pT2(0);
	bool is_neutrino(false);
	bool is_cc(false);
	unordered_map<int, vector<pair<size_t, double>> > primary_pi0_map;
	unordered_map< int, vector<double> > nonprimary_pi0_map;
	double muon_momentum_mag;
	double muon_beam_costheta;
	double pi0_leading_photon_energy;
	TVector3 pi0_leading_photon_dir;
	double pi0_leading_photon_conv_dist;
	double pi0_subleading_photon_energy;
	TVector3 pi0_subleading_photon_dir;
	double pi0_subleading_photon_conv_dist;
	double pi0_costheta;
	double pi0_mass;
	double pi0_momentum_mag;
	double pi0_beam_costheta;

	// Particle loop
	size_t leading_muon_index(0);
	double max_muon_ke(-99999);
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	  const auto & p = obj.particles[i];
	  // Primaries
	  if(p.is_primary)
	  {
	    // Transverse momentum calculation
	    TVector3 _p(p.momentum[0], p.momentum[1], p.momentum[2]);
	    TVector3 pL = _p.Dot(beamdir) * beamdir;
	    TVector3 pT = _p - pL;
	    pT0 += pT[0];
	    pT1 += pT[1];
	    pT2 += pT[2];
		      
	    // Muons
	    if(p.pid == 2)
	    {
	      primary_muon_count++;
	      if(p.ke >= MIN_MUON_ENERGY) primary_muon_count_thresh++;
	      if(p.ke > max_muon_ke)
	      {
		max_muon_ke = p.ke;
		leading_muon_index = i;
	      }
	    }
	    // Pions
	    if(p.pid == 3)
	    {
	      primary_pion_count++;
	      if(p.ke >= MIN_PION_ENERGY) primary_pion_count_thresh++;
	    }
	    // Neutral pions
	    if(p.pdg_code == 22 && p.parent_pdg_code == 111)
	    {
	      primary_pi0_map[p.parent_track_id].push_back({i,p.ke});
	    }
	  } // end primary loop
	  // Nonprimaries
	  else
	  {
	    // Neutral pions
	    if(p.pdg_code == 22 && p.parent_pdg_code == 111)
	    {
	      nonprimary_pi0_map[p.parent_track_id].push_back(p.ke);
	    }
	  } // end nonprimary loop
	} // end particle loop

	// Loop over primary pi0s
	for(auto const & pi0 : primary_pi0_map)
	{
	  int num_primary_photon_daughters(0);
	  int num_primary_photon_daughters_thresh(0);
	  for(auto & pair : pi0.second)
	  {
	    num_primary_photon_daughters++;
	    if(pair.second >= MIN_PHOTON_ENERGY) num_primary_photon_daughters_thresh++;
	  }
	  if(num_primary_photon_daughters == 2) primary_pi0_count++;
	  if(num_primary_photon_daughters_thresh == 2) primary_pi0_count_thresh++;
	}

	// Nonprimary pi0s
	nonprimary_pi0_count = nonprimary_pi0_map.size();

	// Obtain info about signal particles, if they exist
	vector<size_t> pi0_photon_indices;
	if(primary_muon_count_thresh == 1 && primary_pion_count_thresh == 0 && primary_pi0_count_thresh == 1)
	{      
	  // Get leading muon info
	  const auto & muon = obj.particles[leading_muon_index];
	  TVector3 muon_momentum(muon.momentum[0], muon.momentum[1], muon.momentum[2]);
	  double muon_momentum_mag = muon_momentum.Mag();
	  double muon_beam_costheta = muon_momentum.Unit().Dot(beamdir);
	        
	  // Get leading/subleading photon info
	  for(auto const & pi0 : primary_pi0_map)
	  {
	    for(auto pair : pi0.second)
	    {
	      if(pair.second >= MIN_PHOTON_ENERGY) pi0_photon_indices.push_back(pair.first);
	    }
	  }

	  const auto & pi0_photon0 = obj.particles[pi0_photon_indices[0]];
	  const auto & pi0_photon1 = obj.particles[pi0_photon_indices[1]];
	  size_t leading_photon_index;
	  size_t subleading_photon_index;
	  if(pi0_photon0.ke > pi0_photon1.ke)
	  {
	    leading_photon_index = pi0_photon_indices[0];
	    subleading_photon_index = pi0_photon_indices[1];
	  }
	  else
	  {
	    leading_photon_index = pi0_photon_indices[1];
	    subleading_photon_index = pi0_photon_indices[0];
	  }
	  const auto & pi0_leading_photon = obj.particles[leading_photon_index];
	  const auto & pi0_subleading_photon = obj.particles[subleading_photon_index];
	        
	  pi0_leading_photon_energy = pi0_leading_photon.ke;
	  TVector3 pi0_leading_photon_start_point(pi0_leading_photon.start_point[0], pi0_leading_photon.start_point[1], pi0_leading_photon.start_point[2]);
	  TVector3 pi0_leading_photon_dir(pi0_leading_photon.momentum[0], pi0_leading_photon.momentum[1], pi0_leading_photon.momentum[2]);
	  pi0_leading_photon_dir = pi0_leading_photon_dir.Unit();
	  pi0_leading_photon_conv_dist = (vertex - pi0_leading_photon_start_point).Mag();
	  TVector3 pi0_leading_photon_momentum(pi0_leading_photon.momentum[0], pi0_leading_photon.momentum[1], pi0_leading_photon.momentum[2]);

	  pi0_subleading_photon_energy = pi0_subleading_photon.ke;
	  TVector3 pi0_subleading_photon_start_point(pi0_subleading_photon.start_point[0], pi0_subleading_photon.start_point[1], pi0_subleading_photon.start_point[2]);
	  TVector3 pi0_subleading_photon_dir(pi0_subleading_photon.momentum[0], pi0_subleading_photon.momentum[1], pi0_subleading_photon.momentum[2]);
	  pi0_subleading_photon_dir = pi0_subleading_photon_dir.Unit();
	  pi0_subleading_photon_conv_dist = (vertex - pi0_subleading_photon_start_point).Mag();
	  TVector3 pi0_subleading_photon_momentum(pi0_subleading_photon.momentum[0], pi0_subleading_photon.momentum[1],pi0_subleading_photon.momentum[2]);

	  TVector3 pi0_momentum = pi0_leading_photon_momentum + pi0_subleading_photon_momentum;
	  pi0_beam_costheta = pi0_momentum.Unit().Dot(beamdir);

	  pi0_costheta = pi0_leading_photon_dir.Dot(pi0_subleading_photon_dir);
	  pi0_mass = sqrt(2*pi0_leading_photon_energy*pi0_subleading_photon_energy*(1-pi0_costheta));
	        
	  s.muon_momentum_mag = muon_momentum_mag;
	  s.muon_beam_costheta = muon_beam_costheta;
	  s.pi0_leading_photon_energy = pi0_leading_photon_energy;
	  s.pi0_leading_photon_conv_dist = pi0_leading_photon_conv_dist;
	  s.pi0_subleading_photon_energy = pi0_subleading_photon_energy;
	  s.pi0_subleading_photon_conv_dist = pi0_subleading_photon_conv_dist;
	  s.pi0_costheta = pi0_costheta;
	  s.pi0_mass = pi0_mass;
	  s.pi0_momentum_mag = pi0_momentum.Mag();
	  s.pi0_beam_costheta = pi0_beam_costheta;
	} // end signal 
	  
	s.num_primary_muons = primary_muon_count;
	s.num_primary_muons_thresh = primary_muon_count_thresh;
	s.num_primary_pions = primary_pion_count;
	s.num_primary_pions_thresh = primary_pion_count_thresh;
	s.num_primary_pi0s = primary_pi0_count;
	s.num_primary_pi0s_thresh = primary_pi0_count_thresh;
	s.num_nonprimary_pi0s = nonprimary_pi0_count;
	s.transverse_momentum_mag = sqrt(pow(pT0, 2) + pow(pT1, 2) + pow(pT2, 2));
	s.is_fiducial = cuts::fiducial_cut<T>(obj);
	s.has_contained_tracks = cuts::track_containment_cut<T>(obj);
	if(obj.nu_id > -1) is_neutrino = true;
	s.is_neutrino = is_neutrino;
	if(obj.current_type == 0) is_cc = true;
	s.is_cc = is_cc;
	  
	return s;
      }

    /**
     * @brief Structure for reco numu cc pi0 interaction information.
     * @details This structure stores information about reconstructed
     * muon and neutral pion.
     * @tparam T the type of interaction (reco).
     * @param obj the interaction to select on.
     * @return a reco_pi0 structure.
     * @note This structure is intented to be used for the pi02024 analysis.
     */
    template<class T> 
      reco_inter reco_interaction_info(const T & obj)
      {
	// Initialize structure
	reco_inter s;
	  
	// Initialize relevant TVector3s
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
	TVector3 beamdir(0, 0, 1); // BNB

	// Initialize output variables
	double pT0(0), pT1(0), pT2(0);
	size_t leading_muon_index(0);
	size_t leading_photon_index(0);
	size_t subleading_photon_index(0);
	double max_muon_ke(-99999);
	double max_calo_ke0(-99999);
	double max_calo_ke1(-99999);

	// Loop over particles
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	    const auto & p = obj.particles[i];
	    if(!p.is_primary) continue; // Primaries

	    // Transverse momentum calculation
	    TVector3 _p(p.momentum[0], p.momentum[1], p.momentum[2]);
	    TVector3 pL = _p.Dot(beamdir) * beamdir;
	    TVector3 pT = _p - pL;

	    // Muons
	    if(p.pid == 2 && p.csda_ke >= MIN_MUON_ENERGY)
	    {
		if(p.csda_ke > max_muon_ke)
		{
		    max_muon_ke = p.csda_ke;
		    leading_muon_index = i;
		}
	    }
	    // Photons
	    if(p.pid == 0)
	    {
		// Don't use default momentum for photons
		TVector3 _p(p.start_point[0] - vertex[0], p.start_point[1] - vertex[1], p.start_point[2] - vertex[2]);
		_p = p.calo_ke * _p.Unit();
		TVector3 pL = _p.Dot(beamdir) * beamdir;
		TVector3 pT = _p - pL;
		      
		if(p.calo_ke >= MIN_PHOTON_ENERGY)
		{
		    // Leading photon
		    if(p.calo_ke > max_calo_ke0)
		      {
			max_calo_ke0 = p.calo_ke;
			leading_photon_index = i;
		      }
		}
	    }
	    pT0 += pT[0];
	    pT1 += pT[1];
	    pT2 += pT[2];
	} // end particle loop
	  
	// Second particle loop to find subleading photon
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	    const auto & p = obj.particles[i];
	        
            if(p.is_primary && p.pid == 0 && p.calo_ke >= MIN_PHOTON_ENERGY)
	    {
		if(p.calo_ke > max_calo_ke1 && p.calo_ke < max_calo_ke0)
		{
		    max_calo_ke1 = p.calo_ke;
		    subleading_photon_index = i;
		}
	    }
	} // end second particle loop

	// Get leading muon info
	TVector3 muon_momentum;
	muon_momentum.SetX(obj.particles[leading_muon_index].momentum[0]);
	muon_momentum.SetY(obj.particles[leading_muon_index].momentum[1]);
	muon_momentum.SetZ(obj.particles[leading_muon_index].momentum[2]);
	double muon_momentum_mag = muon_momentum.Mag();
	double muon_beam_costheta = muon_momentum.Unit().Dot(beamdir);

	// Get leading photon info
	double pi0_leading_photon_energy = obj.particles[leading_photon_index].calo_ke;
	TVector3 pi0_leading_photon_start_point;
	pi0_leading_photon_start_point.SetX(obj.particles[leading_photon_index].start_point[0]);
	pi0_leading_photon_start_point.SetY(obj.particles[leading_photon_index].start_point[1]);
	pi0_leading_photon_start_point.SetZ(obj.particles[leading_photon_index].start_point[2]);
	double pi0_leading_photon_conv_dist = (vertex - pi0_leading_photon_start_point).Mag();
	TVector3 pi0_leading_photon_dir;
	pi0_leading_photon_dir.SetX(pi0_leading_photon_start_point[0] - vertex[0]);
	pi0_leading_photon_dir.SetY(pi0_leading_photon_start_point[1] - vertex[1]);
	pi0_leading_photon_dir.SetZ(pi0_leading_photon_start_point[2] - vertex[2]);
	pi0_leading_photon_dir = pi0_leading_photon_dir.Unit();
	TVector3 pi0_leading_photon_start_dir;
	pi0_leading_photon_start_dir.SetX(obj.particles[leading_photon_index].start_dir[0]);
	pi0_leading_photon_start_dir.SetY(obj.particles[leading_photon_index].start_dir[1]);
	pi0_leading_photon_start_dir.SetZ(obj.particles[leading_photon_index].start_dir[2]);
	pi0_leading_photon_start_dir = pi0_leading_photon_start_dir.Unit();
	double pi0_leading_photon_cosphi = pi0_leading_photon_dir.Dot(pi0_leading_photon_start_dir);
	TVector3 pi0_leading_photon_momentum = pi0_leading_photon_energy * pi0_leading_photon_dir;

	// Get subleading photon info
	double pi0_subleading_photon_energy = obj.particles[subleading_photon_index].calo_ke;
	TVector3 pi0_subleading_photon_start_point;
	pi0_subleading_photon_start_point.SetX(obj.particles[subleading_photon_index].start_point[0]);
	pi0_subleading_photon_start_point.SetY(obj.particles[subleading_photon_index].start_point[1]);
	pi0_subleading_photon_start_point.SetZ(obj.particles[subleading_photon_index].start_point[2]);
	double pi0_subleading_photon_conv_dist = (vertex - pi0_subleading_photon_start_point).Mag();
	TVector3 pi0_subleading_photon_dir;
	pi0_subleading_photon_dir.SetX(pi0_subleading_photon_start_point[0] - vertex[0]);
	pi0_subleading_photon_dir.SetY(pi0_subleading_photon_start_point[1] - vertex[1]);
	pi0_subleading_photon_dir.SetZ(pi0_subleading_photon_start_point[2] - vertex[2]);
	pi0_subleading_photon_dir = pi0_subleading_photon_dir.Unit();
	TVector3 pi0_subleading_photon_start_dir;
	pi0_subleading_photon_start_dir.SetX(obj.particles[subleading_photon_index].start_dir[0]);
	pi0_subleading_photon_start_dir.SetY(obj.particles[subleading_photon_index].start_dir[1]);
	pi0_subleading_photon_start_dir.SetZ(obj.particles[subleading_photon_index].start_dir[2]);
	pi0_subleading_photon_start_dir = pi0_subleading_photon_start_dir.Unit();
	double pi0_subleading_photon_cosphi = pi0_subleading_photon_dir.Dot(pi0_subleading_photon_start_dir);
	TVector3 pi0_subleading_photon_momentum = pi0_subleading_photon_energy * pi0_subleading_photon_dir;

	// Get neutral pion info
	double pi0_cos_opening_angle = pi0_leading_photon_dir.Dot(pi0_subleading_photon_dir);
	double pi0_mass = sqrt(2*pi0_leading_photon_energy*pi0_subleading_photon_energy*(1-pi0_cos_opening_angle));
	TVector3 pi0_momentum = pi0_leading_photon_momentum + pi0_subleading_photon_momentum;
	double pi0_momentum_mag = pi0_momentum.Mag();
	double pi0_beam_costheta = pi0_momentum.Unit().Dot(beamdir);

	// Fill struct
	s.transverse_momentum_mag = sqrt(pow(pT0, 2) + pow(pT1, 2) + pow(pT2, 2));
	s.muon_momentum_mag = muon_momentum_mag;
	s.muon_beam_costheta = muon_beam_costheta;
	s.pi0_leading_photon_energy = pi0_leading_photon_energy;
	s.pi0_leading_photon_cosphi = pi0_leading_photon_cosphi;
	s.pi0_leading_photon_conv_dist = pi0_leading_photon_conv_dist;
	s.pi0_subleading_photon_energy = pi0_subleading_photon_energy;
	s.pi0_subleading_photon_cosphi = pi0_subleading_photon_cosphi;
	s.pi0_subleading_photon_conv_dist = pi0_subleading_photon_conv_dist;
	s.pi0_costheta = pi0_cos_opening_angle;
	s.pi0_mass = pi0_mass;
	s.pi0_momentum_mag = pi0_momentum_mag;
	s.pi0_beam_costheta = pi0_beam_costheta;

	return s;
      }

}
#endif // UTILITIES_PI02024_H
