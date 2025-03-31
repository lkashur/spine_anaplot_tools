/**
 * @file utilities_pi0ana_nc.h
 * @brief Header file for definitions of utility functions for supporting
 * analysis variables and cuts.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author lkashur@colostate.edu
 */
#ifndef UTILITIES_PI0ANA_NC_H
#define UTILITIES_PI0ANA_NC_H

#include <iostream>
#include <vector>
#include <TVector3.h>
#include "include/cuts.h"

#define MIN_MUON_ENERGY 25
#define MIN_PROTON_ENERGY 40
#define MIN_PHOTON_ENERGY 40
#define MIN_OTHER_ENERGY 25

struct truth_inter_nc {
  int num_primary_muons;
  int num_primary_muons_thresh;
  int num_primary_pions;
  int num_primary_pions_thresh;
  int num_primary_pi0s;
  int num_primary_pi0s_thresh;
  int num_nonprimary_pi0s;
  int num_primary_protons;
  int num_primary_protons_thresh;
  double transverse_momentum_mag;
  bool is_fiducial;
  bool has_contained_tracks;
  bool is_neutrino;
  bool is_cc;
  double pi0_leading_photon_energy;
  double pi0_leading_photon_conv_dist;
  double pi0_subleading_photon_energy;
  double pi0_subleading_photon_conv_dist;
  double pi0_photons_costheta;
  double pi0_mass;
  double pi0_momentum_mag;
  double pi0_beam_costheta;
};

struct reco_inter_nc {
  double transverse_momentum_mag;
  double pi0_leading_photon_energy;
  double pi0_leading_photon_cosphi;
  double pi0_leading_photon_conv_dist;
  double pi0_subleading_photon_energy;
  double pi0_subleading_photon_cosphi;
  double pi0_subleading_photon_conv_dist;
  double pi0_photons_costheta;
  double pi0_mass;
  double pi0_momentum_mag;
  double pi0_beam_costheta;
};



/**
 * @namespace utilities_pi0ana_nc
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
namespace utilities_pi0ana_nc
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
	    if(PIDFUNC(p) == 0 && p.ke >= MIN_PHOTON_ENERGY) passes = true; // Photons
	    if(PIDFUNC(p) == 1 && p.ke >= MIN_OTHER_ENERGY) passes = true; // Electrons
	    if(PIDFUNC(p) == 2 && p.ke >= MIN_MUON_ENERGY) passes = true; // Muons
	    if(PIDFUNC(p) == 3 && p.ke >= MIN_OTHER_ENERGY) passes = true; // Pions
	    if(PIDFUNC(p) == 4 && p.ke >= MIN_PROTON_ENERGY) passes = true; // Protons
	    if(PIDFUNC(p) == 5 && p.ke >= MIN_OTHER_ENERGY) passes = true; // Kaons
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
              ++counts[PIDFUNC(p)];
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
     * @note This structure is intented to be used for the pi0ana analysis. 
     */
    template<class T> 
      truth_inter_nc truth_interaction_info(const T & obj)
      {
	// Initialize struct
	truth_inter_nc s;
	  
	// Initialize relevant TVector3s
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
	TVector3 beamdir;
        if constexpr(!BEAM_IS_NUMI){
            beamdir.SetX(0);
            beamdir.SetY(0);
            beamdir.SetZ(1);
          }
        else{
          beamdir.SetX(315.120380 + vertex[0]);
          beamdir.SetY(33.644912 + vertex[1]);
          beamdir.SetZ(733.632532 + vertex[2]);
        }
	//TVector3 beamdir(BEAMDIR);
	  
	// Initialize output variables
	int primary_muon_count(0);
	int primary_muon_count_thresh(0);
	int primary_pion_count(0);
	int primary_pion_count_thresh(0);
	int primary_pi0_count(0);
	int primary_pi0_count_thresh(0);
	int nonprimary_pi0_count(0);
	int primary_proton_count(0);
	int primary_proton_count_thresh(0);
	double pT0(0), pT1(0), pT2(0);
	bool is_neutrino(false);
	bool is_cc(false);
	unordered_map<int, vector<pair<size_t, double>> > primary_pi0_map;
	unordered_map<int, vector<pair<size_t, double>> > nonprimary_pi0_map;
	double muon_momentum_mag;
	double muon_csda_ke;
	double muon_calo_ke;
	double muon_beam_costheta;
	double pi0_leading_photon_energy;
	TVector3 pi0_leading_photon_dir;
	double pi0_leading_photon_conv_dist;
	double pi0_subleading_photon_energy;
	TVector3 pi0_subleading_photon_dir;
	double pi0_subleading_photon_conv_dist;
	double pi0_photons_costheta;
	double pi0_mass;
	TVector3 pi0_momentum;
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
	    if(PIDFUNC(p) == 2)
	    {
	      primary_muon_count++;
	      if(p.ke >= MIN_MUON_ENERGY) primary_muon_count_thresh++;
	      //std::cout << "Muon KE: " << p.ke << std::endl;
	      if(p.ke > max_muon_ke)
	      {
		max_muon_ke = p.ke;
		leading_muon_index = i;
	      }
	    }
	    // Pions
	    if(PIDFUNC(p) == 3)
	    {
	      primary_pion_count++;
	      if(p.ke >= MIN_OTHER_ENERGY) primary_pion_count_thresh++;
	    }

	    // Protons
	    if(PIDFUNC(p) == 4)
	    {
	      primary_proton_count++;
	      if(p.ke >= MIN_PROTON_ENERGY) primary_proton_count_thresh++;
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
	      nonprimary_pi0_map[p.parent_track_id].push_back({i,-5});
	    }
	  } // end nonprimary loop
	} // end particle loop

	// Loop over primary pi0s
	vector<int> bad_primary_pi0_ids;
	for(auto const & pi0 : primary_pi0_map)
	  {
	    // Loop over daughters of each pi0
	    int num_primary_photon_daughters(0);
	    for(auto & daughter : pi0.second)
	      {
		num_primary_photon_daughters++;
	      }
	    if(num_primary_photon_daughters != 2) bad_primary_pi0_ids.push_back(pi0.first);
	  }

	// Remove dalitz and single photon pi0s from map
	for(size_t i=0; i<bad_primary_pi0_ids.size(); i++)
	{
	    primary_pi0_map.erase(bad_primary_pi0_ids[i]);
	}
        primary_pi0_count = primary_pi0_map.size();

	// Primary pi0s above threshold
	vector<int> subthresh_primary_pi0_ids;
	for(auto const & pi0 : primary_pi0_map)
	  {
	    int num_photon_daughters_thresh(0);
	    for(auto & daughter : pi0.second)
	      {
		if(daughter.second >= MIN_PHOTON_ENERGY) num_photon_daughters_thresh++;
	      }
	    
	    if(num_photon_daughters_thresh == 2)
	      {
		primary_pi0_count_thresh++;
	      }
	    else
	      {
		subthresh_primary_pi0_ids.push_back(pi0.first);
	      }
	  }

	// Remove subthreshold pi0s from our map
	for(size_t i=0; i<subthresh_primary_pi0_ids.size(); i++)
	  {
            primary_pi0_map.erase(subthresh_primary_pi0_ids[i]);
	  }

	// Nonprimary pi0s
	nonprimary_pi0_count = nonprimary_pi0_map.size();

	// Obtain info about signal particles, if they exist
	if(primary_muon_count_thresh == 0 && primary_pion_count_thresh == 0 && primary_pi0_count_thresh == 1)
	{      
	        
	  // Get leading/subleading photon info
	  vector<size_t> pi0_photon_indices;
	  for(auto const & pi0 : primary_pi0_map)
	  {
	    for(auto daughter : pi0.second)
	    {
	      pi0_photon_indices.push_back(daughter.first);
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
	        
	  pi0_leading_photon_energy = pi0_leading_photon.ke; //.ke
	  TVector3 pi0_leading_photon_start_point(pi0_leading_photon.start_point[0], pi0_leading_photon.start_point[1], pi0_leading_photon.start_point[2]);
	  TVector3 pi0_leading_photon_dir(pi0_leading_photon.momentum[0], pi0_leading_photon.momentum[1], pi0_leading_photon.momentum[2]);
	  pi0_leading_photon_dir = pi0_leading_photon_dir.Unit();
	  pi0_leading_photon_conv_dist = (vertex - pi0_leading_photon_start_point).Mag();
	  TVector3 pi0_leading_photon_momentum(pi0_leading_photon.momentum[0], pi0_leading_photon.momentum[1], pi0_leading_photon.momentum[2]);

	  pi0_subleading_photon_energy = pi0_subleading_photon.ke; // .ke
	  TVector3 pi0_subleading_photon_start_point(pi0_subleading_photon.start_point[0], pi0_subleading_photon.start_point[1], pi0_subleading_photon.start_point[2]);
	  TVector3 pi0_subleading_photon_dir(pi0_subleading_photon.momentum[0], pi0_subleading_photon.momentum[1], pi0_subleading_photon.momentum[2]);
	  pi0_subleading_photon_dir = pi0_subleading_photon_dir.Unit();
	  pi0_subleading_photon_conv_dist = (vertex - pi0_subleading_photon_start_point).Mag();
	  TVector3 pi0_subleading_photon_momentum(pi0_subleading_photon.momentum[0], pi0_subleading_photon.momentum[1], pi0_subleading_photon.momentum[2]);

	  pi0_momentum = pi0_leading_photon_momentum + pi0_subleading_photon_momentum;
	  pi0_beam_costheta = pi0_momentum.Unit().Dot(beamdir);

	  pi0_photons_costheta = pi0_leading_photon_dir.Dot(pi0_subleading_photon_dir);
	  pi0_mass = sqrt(2*pi0_leading_photon_energy*pi0_subleading_photon_energy*(1-pi0_photons_costheta));
	  
	  s.pi0_leading_photon_energy = pi0_leading_photon_energy;
	  s.pi0_leading_photon_conv_dist = pi0_leading_photon_conv_dist;
	  s.pi0_subleading_photon_energy = pi0_subleading_photon_energy;
	  s.pi0_subleading_photon_conv_dist = pi0_subleading_photon_conv_dist;
	  s.pi0_photons_costheta = pi0_photons_costheta;
	  s.pi0_mass = pi0_mass;
	  s.pi0_momentum_mag = pi0_momentum.Mag();
	  s.pi0_beam_costheta = pi0_beam_costheta;
	} // end signal 
	else
	{
          s.pi0_leading_photon_energy = -5;
          s.pi0_leading_photon_conv_dist = -5;
          s.pi0_subleading_photon_energy = -5;
          s.pi0_subleading_photon_conv_dist = -5;
          s.pi0_photons_costheta = -5;
          s.pi0_mass = -5;
          s.pi0_momentum_mag = -5;
          s.pi0_beam_costheta = -5;
	}
	
	s.num_primary_muons = primary_muon_count;
	s.num_primary_muons_thresh = primary_muon_count_thresh;
	s.num_primary_pions = primary_pion_count;
	s.num_primary_pions_thresh = primary_pion_count_thresh;
	s.num_primary_pi0s = primary_pi0_count;
	s.num_primary_pi0s_thresh = primary_pi0_count_thresh;
	s.num_nonprimary_pi0s = nonprimary_pi0_count;
	s.num_primary_protons = primary_proton_count;
	s.num_primary_protons_thresh = primary_proton_count_thresh;
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
     * @note This structure is intented to be used for the pi0ana analysis.
     */
    template<class T> 
      reco_inter_nc reco_interaction_info(const T & obj)
      {
	// Initialize structure
	reco_inter_nc s;
	  
	// Initialize relevant TVector3s
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
        TVector3 beamdir;
        if constexpr(!BEAM_IS_NUMI){
            beamdir.SetX(0);
            beamdir.SetY(0);
            beamdir.SetZ(1);
          }
        else{
          beamdir.SetX(315.120380 + vertex[0]);
          beamdir.SetY(33.644912 + vertex[1]);
          beamdir.SetZ(733.632532 + vertex[2]);
        }
	//TVector3 beamdir(BEAMDIR);

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

	    // Photons
	    if(PIDFUNC(p) == 0)
	    {
		// Don't use default momentum for photons
		TVector3 _p(p.start_point[0] - vertex[0], p.start_point[1] - vertex[1], p.start_point[2] - vertex[2]);
		_p = p.calo_ke * _p.Unit();
		TVector3 pL = _p.Dot(beamdir) * beamdir;
		TVector3 pT = _p - pL;
		 
		// Leading photon
		if(p.calo_ke > max_calo_ke0)
		{
		    max_calo_ke0 = p.calo_ke;
		    leading_photon_index = i;
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
	        
            if(p.is_primary && PIDFUNC(p) == 0 && p.calo_ke >= MIN_PHOTON_ENERGY)
	    {
		if(p.calo_ke > max_calo_ke1 && p.calo_ke < max_calo_ke0)
		{
		    max_calo_ke1 = p.calo_ke;
		    subleading_photon_index = i;
		}
	    }
	} // end second particle loop

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
	double pi0_photons_costheta = pi0_leading_photon_dir.Dot(pi0_subleading_photon_dir);
	double pi0_mass = sqrt(2*pi0_leading_photon_energy*pi0_subleading_photon_energy*(1-pi0_photons_costheta));
	TVector3 pi0_momentum = pi0_leading_photon_momentum + pi0_subleading_photon_momentum;
	double pi0_momentum_mag = pi0_momentum.Mag();
	double pi0_beam_costheta = pi0_momentum.Unit().Dot(beamdir);

	// Fill struct
	s.transverse_momentum_mag = sqrt(pow(pT0, 2) + pow(pT1, 2) + pow(pT2, 2));
	s.pi0_leading_photon_energy = pi0_leading_photon_energy;
	s.pi0_leading_photon_cosphi = pi0_leading_photon_cosphi;
	s.pi0_leading_photon_conv_dist = pi0_leading_photon_conv_dist;
	s.pi0_subleading_photon_energy = pi0_subleading_photon_energy;
	s.pi0_subleading_photon_cosphi = pi0_subleading_photon_cosphi;
	s.pi0_subleading_photon_conv_dist = pi0_subleading_photon_conv_dist;
	s.pi0_photons_costheta = pi0_photons_costheta;
	s.pi0_mass = pi0_mass;
	s.pi0_momentum_mag = pi0_momentum_mag;
	s.pi0_beam_costheta = pi0_beam_costheta;

	return s;
      }

}
#endif // UTILITIES_PI0ANA_NC_H
