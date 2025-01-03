/**
 * @file vars_pi02024.h
 * @brief Header file for definitions of analysis variables specific to the
 * pi02024 analysis.
 * @details This file contains definitions of analysis variables which can be
 65;6003;1c* used to extract information from interactions specific to the pi02024
 * analysis. Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author lkashur@colostate.edu
 */
#ifndef VARS_PI02024_H
#define VARS_PI02024_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/utilities.h"
#include "include/cuts.h"
#include "include/pi02024/cuts_pi02024.h"

/**
 * @namespace vars::pi02024
 * @brief Namespace for organizing variables specific to the pi02024 analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the pi02024 analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::pi02024
{
    /**
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the 
     * three categories.
     * 0: Signal (contained and fiducial)
     * 1: Signal (not contained or not fiducial)
     * 2: Other nu
     * 3: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the interaction.
     */
    double category(const caf::SRInteractionTruthDLPProxy & obj)
    {
      // Cosmic background
      double cat(3);

      // Signal   
      if(cuts::pi02024::signal_1mu0pi1pi0(obj))
      {
	if(cuts::fiducial_cut(obj) && cuts::track_containment_cut(obj))
	{
	  cat = 0;
	}
	else cat = 1;
      }
      // Neutrino Background              
      else if(cuts::pi02024::other_nu_1mu0pi1pi0(obj))
      {
	cat = 2;
      }
      return cat;
    }

    /**
     * @brief Variable for enumerating cut type.
     * @details This variable provides a basic categorization of cuts
     * using only signal and sideband as the two cateogories.
     * 1: Signal
     * 2: Sideband
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the cut. 
     */
    template<class T>
        double cut_type(const T & obj)
        {
	  // Signal
	  double cat(1);

	  return cat;
        }

    /**
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 1: Signal
     * 2: Signal (OOPS)
     * 3: Other nu
     * 4: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the interaction.
     */
    double is_signal(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter s = utilities_pi02024::truth_interaction_info(obj);
	
	// Cosmic
	uint16_t cat(4);
	
	// Nu
	if(s.is_neutrino)
	{
	    // Signal
	    if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial && s.has_contained_tracks) cat = 1;

	    // Signal (OOPS)
	    else if( (s.num_primary_muons == 1 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.is_cc) && (s.num_primary_muons_thresh != 1 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1 || !s.is_fiducial || !s.has_contained_tracks) ) cat = 2;
	    
	    // Other nu
	    else cat = 3;
	}
	
	return cat;
    }

    /**
     * @brief Variable for enumerating interaction topological categories.
     * @details This variable provides a basic categorization of interaction
     * topologies using the following cateogories:
     * 0: 1mu0pi1pi0 (signal)
     * 1: 1mu0pi1pi0 (OOPs)
     * 2: 1muNpi1pi0
     * 3: 1muNpi0pi0
     * 4: 1muNpi0
     * 5: NC 1pi0
     * 6: Other nu
     * 7: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topological category of the interaction.
     */
    double category_topology(const caf::SRInteractionTruthDLPProxy & obj)
    {
      truth_inter s = utilities_pi02024::truth_interaction_info(obj);

      // Cosmic                                                               
      uint16_t cat(7);

      // Neutrino                                                                                              
      if(s.is_neutrino)
      {
	// 1mu0pi1pi0 (in-phase)
	if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial && s.has_contained_tracks) cat = 0;
	// 1mu0pi1pi0 (out-of-phase)
	else if( (s.num_primary_muons == 1 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.is_cc) && (s.num_primary_muons_thresh != 1 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1 || !s.is_fiducial || !s.has_contained_tracks) ) cat = 1;
	// 1muNpi1pi0
	else if(s.num_primary_muons == 1 && s.num_primary_pions > 0 && s.num_primary_pi0s == 1 && s.is_cc) cat = 2;
	// 1muNpi0pi0
	else if(s.num_primary_muons == 1 && s.num_primary_pions > 0 && s.num_primary_pi0s == 0 && s.is_cc) cat = 3;
	// 1muNpi0
	else if(s.num_primary_muons == 1 && s.num_primary_pi0s > 1 && s.is_cc) cat = 4;
	// NC 1pi0
	else if(s.num_primary_muons == 0 && s.num_primary_pi0s == 1 && !s.is_cc) cat = 5;
	// Other nu
	else cat = 6;
      }
      return cat;
    }

    /**
     * @brief Variable for leading muon momentum magnitude.
     * @details Variable for momentum of the leading muon
     * candidate [MeV/c].
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the muon momentum magnitude.
     */
    template<class T> 
        double muon_momentum_mag(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_pi02024::truth_interaction_info(obj);
			       return s.muon_momentum_mag;
			   }
	    else
            {
	        reco_inter s = utilities_pi02024::reco_interaction_info(obj);
		return s.muon_momentum_mag;
            }
        }

    /**
     * @brief Variable for leading muon angle with beam.
     * @details Variable for the cosine of the angle between
     * the interaction's leading muon and the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the leading muon's angle w.r.t. beam.
     */
    template<class T>
        double muon_beam_costheta(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                           {
			       truth_inter s = utilities_pi02024::truth_interaction_info(obj);
			       return s.muon_beam_costheta;
                           }
	    else
	    {
	      reco_inter s = utilities_pi02024::reco_interaction_info(obj);
	      return s.muon_beam_costheta;
	    }
        }

    /**
     * @brief Variable for pi0 leading photon energy.
     * @details Variable for pi0 leading photon energy
     * [MeV], as calculated with p.calo_ke attribute.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 leading photon energy.
     */
    template<class T>
        double pi0_leading_photon_energy(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_pi02024::truth_interaction_info(obj);
			       return s.pi0_leading_photon_energy;
			   }
            else
	    {
		reco_inter s = utilities_pi02024::reco_interaction_info(obj);
		return s.pi0_leading_photon_energy;
	    }
	}

    /**
     * @brief Variable for pi0 leading photon conversion distance.
     * @details Variable for pi0 leading photon conversion distance
     * [cm], as calculated using interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 leading photon conversion distance.
     */
    template<class T>
      double pi0_leading_photon_conv_dist(const T & obj)
      {
	  if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			 {
			     truth_inter s = utilities_pi02024::truth_interaction_info(obj);
			     return s.pi0_leading_photon_conv_dist;
			 }
	  else
	  {
	      reco_inter s = utilities_pi02024::reco_interaction_info(obj);
	      return s.pi0_leading_photon_conv_dist;
	  }
      }

    /**
     * @brief Variable for pi0 subleading photon energy.
     * @details Variable for pi0 subleading photon energy
     * [MeV], as calculated with p.calo_ke attribute.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 subleading photon energy. 
     */
    template<class T>
        double pi0_subleading_photon_energy(const T & obj)
	{
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                             {
			       truth_inter s = utilities_pi02024::truth_interaction_info(obj);
			       return s.pi0_subleading_photon_energy;
                             }
            else
	    {
	        reco_inter s = utilities_pi02024::reco_interaction_info(obj);
                return s.pi0_subleading_photon_energy;
	    }
	}

    /**
     * @brief Variable for pi0 subleading photon conversion distance.
     * @details Variable for pi0 subleading photon conversion distance
     * [cm], as calculated using interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco). 
     * @param obj the interaction to apply the variable on.
     * @return the pi0 subleading photon conversion distance. 
     */
    template<class T>
      double pi0_subleading_photon_conv_dist(const T & obj)
      {
	  if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		         {
			     truth_inter s = utilities_pi02024::truth_interaction_info(obj);
			     return s.pi0_subleading_photon_conv_dist;
		         }
	  else
          {
	      reco_inter s = utilities_pi02024::reco_interaction_info(obj);
	      return s.pi0_subleading_photon_conv_dist;
          }
      }

    /**
     * @brief Variable for neutral pion mass.
     * @details Variable for neutral pion mass [MeV/c^2], as calculated
     * with photon energies and opening angle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion mass.
     */
    template<class T>
        double pi0_mass(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                           {
			       truth_inter s = utilities_pi02024::truth_interaction_info(obj);
                               return s.pi0_mass;
                           }
	    else
	    {
	        reco_inter s = utilities_pi02024::reco_interaction_info(obj);
	        return s.pi0_mass;
	    }
	} 
}
#endif // VARS_PI02024_H
