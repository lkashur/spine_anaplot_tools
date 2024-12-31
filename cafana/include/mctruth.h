/**
 * @file mctruth.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 */
#ifndef MCTRUTH_H
#define MCTRUTH_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing variables which act on true interactions.
 * @details This namespace is intended to be used for organizing variables
 * which act on true interactions. Each variable is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a double.
 */
namespace mctruth
{
    /**
     * @brief Variable for the true neutrino energy.
     * @details This variable is intended to provide the true energy of the
     * parent neutrino that produced the interaction.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy.
     */
    double true_neutrino_energy(const caf::Proxy<caf::SRTrueInteraction> & obj) { return obj.E; }

    /**
     * @brief Variable for the true neutrino baseline.
     * @details This variable is intended to provide the true baseline of the
     * parent neutrino that produced the interaction.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino baseline.
     */
    double true_neutrino_baseline(const caf::Proxy<caf::SRTrueInteraction> & obj) { return obj.baseline; }

    /**
     * @brief Variable for the true neutrino PDG code.
     * @details This variable is intended to provide the true PDG code of the
     * parent neutrino that produced the interaction.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    double true_neutrino_pdg(const caf::Proxy<caf::SRTrueInteraction> & obj) { return obj.pdg; }

    /**
     * @brief Variable for the true neutrino current value.
     * @details This variable is intended to provide the true current value of
     * the parent neutrino that produced the interaction.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino current value.
     */
    double true_neutrino_cc(const caf::Proxy<caf::SRTrueInteraction> & obj) { return obj.iscc; }

    /**
     * @brief Variable for the interaction mode of the interaction.
     * @details This variable is intended to provide the interaction mode of the
     * interaction. This is based on the GENIE interaction mode enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction mode.
     */
    double interaction_mode(const caf::Proxy<caf::SRTrueInteraction> & obj) { return obj.genie_mode; }

    /**
     * @brief Variable for the interaction type of the interaction.
     * @details This variable is intended to provide the interaction type of the
     * interaction. This is based on the GENIE interaction type enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction type.
     */
    double interaction_type(const caf::Proxy<caf::SRTrueInteraction> & obj) { return obj.genie_inttype; }
} // namespace mctruth
#endif