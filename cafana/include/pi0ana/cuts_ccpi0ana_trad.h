/**
 * @file cuts_ccpi0ana_trad.h
 * @brief Header file for definitions of analysis cuts specific to the muonana
 * analysis.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the ccpi0ana analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author lkashur@colostate.edu
*/
#ifndef CUTS_CCPI0ANA_TRAD_H
#define CUTS_CCPI0ANA_TRAD_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "utilities_ccpi0ana_trad.h"

/**
 * @namespace cuts::ccpi0ana_trad
 * @brief Namespace for organizing cuts specific to the ccpi0ana analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the ccpi0ana analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::ccpi0ana_trad
{
    /**
     * @brief Apply data cut.
     * @details This cut asserts the interaction is from data (not MC).
     * @param sr the Standard Record.
     * @return true if interaction is from data.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
  //bool is_data(const caf::Proxy<caf::StandardRecord> sr) 
  //{
  //bool _is_data(sr.ndlp_true == 0);
  //return _is_data;
  //}
  
    /**
     * @brief Apply a 1mu0pi2gamma topological (final state) cut.
     * @details The interaction must have a topology matching 1mu0pi2gamma as defined by
     * the conditions in the @ref count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1mu0pi2gamma topology.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    template<class T>
        bool topological_1mu0pi2gamma_cut(const T & obj)
        {
	    std::vector<uint32_t> c(utilities_ccpi0ana_trad::count_primaries(obj));
            return c[0] == 2 && c[2] == 1 && c[3] == 0;
        }

    template<class T>
      bool zero_charged_pions_cut(const T & obj)
      {
	std::vector<uint32_t> c(utilities_ccpi0ana_trad::count_primaries(obj));
	return c[3] == 0;
      }

    template<class T>
      bool one_muon_cut(const T & obj)
      {
	std::vector<uint32_t> c(utilities_ccpi0ana_trad::count_primaries(obj));
	return c[2] == 1;
      }

    template<class T>
      bool two_photons_cut(const T & obj)
      {
	std::vector<uint32_t> c(utilities_ccpi0ana_trad::count_primaries(obj));
        return c[0] == 2;
      }

    template<class T>
      bool two_or_three_photons_cut(const T & obj)
      {
	std::vector<uint32_t> c(utilities_ccpi0ana_trad::count_primaries(obj));
	return c[0] > 1 & c[0] < 4;
      }

      
    /**
     * @brief Apply pi0 mass cut.
     * @details This function applies a cut on the invariant diphoton mass
     * of the interactions two most energetic photons.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the pi0 mass cut.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    template<class T>
        bool pi0_mass_cut(const T & obj)
        {
	  reco_inter_trad s = utilities_ccpi0ana_trad::reco_interaction_info(obj);
	  return s.pi0_mass < 400;
	}

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), 1mu0pi2gamma
     * topological, and pi0 mass cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), 1mu0pi2gamma topological, and pi0 mass cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, 1mu0pi2gamma topological, and pi0 mass cut.
     * @note This cut is intended to be used for the ccpi0ana analysis. 
     */
    template<class T>
        bool all_1mu0pi2gamma_cut(const T & obj) {return fiducial_cut<T>(obj) && track_containment_cut<T>(obj) && flash_cut<T>(obj) && topological_1mu0pi2gamma_cut<T>(obj) && pi0_mass_cut<T>(obj);}

    /**
     * @brief Apply a cut to select the 1mu0pi1pi0 signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction.  This is the "true" 1mu0pi1pi0 signal.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu0pi1pi0 topological cut.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    bool signal_1mu0pi1pi0(const caf::SRInteractionTruthDLPProxy & obj)
        {
	  truth_inter_trad s = utilities_ccpi0ana_trad::truth_interaction_info(obj);
	  return s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_neutrino;
        }

    /**
     * @brief Apply a cut to select other (non 1mu0pi1pi0) neutrino interactions.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction.  This is the "true" non-signal cut.
     * (neutrino inteaction, but not 1mu0pi1pi0)
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment, 
     * and non-1mu0pi1pi0 topological cut.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    bool other_nu_1mu0pi1pi0(const caf::SRInteractionTruthDLPProxy & obj)
        {
	  truth_inter_trad s = utilities_ccpi0ana_trad::truth_interaction_info(obj);
	  return !(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc) && s.is_neutrino;
        }

}
#endif // CUTS_CCPI0ANA_TRAD_H
