/**
 * @file vars_pi0ana_phase.h
 * @brief Header file for definitions of analysis variables specific to the
 * pi0ana analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions specific to the pi0ana
 * analysis. Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author lkashur@colostate.edu
 */
#ifndef VARS_PI0ANA_PHASE_H
#define VARS_PI0ANA_PHASE_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include <iostream>
#include "include/utilities.h"
#include "include/cuts.h"
#include "include/pi0ana/cuts_pi0ana_phase.h"

/**
 * @namespace vars::pi0ana_phase
 * @brief Namespace for organizing variables specific to the pi0ana analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the pi0ana analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::pi0ana_phase
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
      if(cuts::pi0ana_phase::signal_1mu0pi1pi0(obj))
      {
	if(cuts::fiducial_cut(obj))
	{
	  cat = 0;
	}
	else cat = 1;
      }
      // Neutrino Background              
      else if(cuts::pi0ana_phase::other_nu_1mu0pi1pi0(obj))
      {
	cat = 2;
      }
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
      truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);

      // Cosmic                                                                                                              
      uint16_t cat(4);

      // Nu                                                                                                                         
      if(s.is_neutrino)
        {
	  // Signal
	  if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 1;

	  // Signal (OOPS)
	  else if( (s.num_primary_muons == 1 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.is_cc) && (s.num_primary_muons_thresh != 1 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1 || !s.is_fiducial) ) cat = 2;

	  // Other nu
	  else cat = 3;
        }

      return cat;
    }

    double category_topology_simple(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);

	// Cosmic
	uint16_t cat(5);
	
	// Neutrino
	if(s.is_neutrino)
	{
	    // 1mu 0pi 1pi0 (in-phase)
	    if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 0;
	    // 1mu 0pi Npi0
	    else if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh > 1 && s.is_cc) cat = 1;
	    // 1mu Npi
	    else if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh > 0 && s.is_cc) cat = 2;
	    // Other CC nu
	    else if(s.is_cc) cat = 3;
	    // NC nu
	    else cat = 4;
	}
	return cat;
    }


    double category_topology(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);

	// Cosmic
	uint16_t cat(7);

	// Neutrino
	if(s.is_neutrino)
	{
	  // 1mu0pi1pi0 (in-phase)
	  if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 0;
	  // 1mu0pi1pi0 (OOPS)
	  else if( (s.num_primary_muons == 1 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.is_cc) && (s.num_primary_muons_thresh != 1 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1 || !s.is_fiducial) ) cat = 1;
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

    template<class T>
        double is_not_data(const T & obj)
        {
	    double cat(0);
	    return cat;
	}
 
    template<class T>
        double is_data(const T & obj)
	{
            double cat(1);
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
			       truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			       return s.muon_momentum_mag;
			   }
	    else
            {
	        reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
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
			       truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			       return s.muon_beam_costheta;
                           }
	    else
	    {
	      reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
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
			       truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			       return s.pi0_leading_photon_energy;
			   }
            else
	    {
		reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
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
			     truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			     return s.pi0_leading_photon_conv_dist;
			 }
	  else
	  {
	      reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
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
			       truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			       return s.pi0_subleading_photon_energy;
                             }
            else
	    {
	        reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
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
			     truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			     return s.pi0_subleading_photon_conv_dist;
		         }
	  else
          {
	      reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
	      return s.pi0_subleading_photon_conv_dist;
          }
      }

    /**
     * @brief Variable for neutral pion momentum magnitude.
     * @details Variable for momentum of the neutral pion
     * candidate [MeV/c].
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion momentum.
     */
    template<class T>
        double pi0_momentum_mag(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			   return s.pi0_momentum_mag;
		       }
	    else
	    {
	        reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
		return s.pi0_momentum_mag;
	    } 
        }

    /**
     * @brief Variable for neutral pion angle with beam.
     * @details Variable for the cosine of the angle between
     * the interaction's neutral pion and the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion angle w.r.t. beam.
     */
    template<class T>
        double pi0_beam_costheta(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			   return s.pi0_beam_costheta;
		       }
	    else
	    {
	        reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
		return s.pi0_beam_costheta;
	    }
        }

    /**
     * @brief Variable for neutral pion (photons) opening angle.
     * @details Variable for the openign angle between neutral pion
     * photons, as calculated using the interaction vertex and shower
     * start points.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     */
    template<class T>
        double pi0_photons_costheta(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
			       return s.pi0_photons_costheta;
			   }
	    else
	    {
	        reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
		return s.pi0_photons_costheta;
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
			       truth_inter_phase s = utilities_pi0ana_phase::truth_interaction_info(obj);
                               return s.pi0_mass;
                           }
	    else
	    {
	        reco_inter_phase s = utilities_pi0ana_phase::reco_interaction_info(obj);
	        return s.pi0_mass;
	    }
	} 

    /**
     * @brief Variable for total visible energy of interaction.        
     * @details This function calculates the total visible energy of the   
     * interaction by summing the energy of all particles that are identified                  
     * as counting towards the final state of the interaction.                        
     * @tparam T the type of interaction (true or reco).                         
     * @param obj interaction to apply the variable on.                     
     * @return the total visible energy of the interaction.                                                                         
     */
    template<class T>
        double visible_energy(const T & obj)
        {
	    double energy(0);
	    for(const auto & p : obj.particles)
	    {
	        if(utilities_pi0ana_phase::final_state_signal(p))
		{
		    energy += pvars::energy(p);
		    if(PIDFUNC(p) == 4) energy -= pvars::mass(p) - PROTON_BINDING_ENERGY;
		}
	    }
	    return energy/1000.0;
        }

}
#endif // VARS_PI0ANA_PHASE_H
