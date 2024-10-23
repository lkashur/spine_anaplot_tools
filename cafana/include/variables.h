/**
 * @file variables.h
 * @brief Header file for definitions of analysis variables.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions. Each variable is implemented
 * as a function which takes an interaction object as an argument and returns a
 * double. These are the building blocks for producing high-level plots of the
 * selected interactions.
 * @author mueller@fnal.gov
*/
#ifndef VARIABLES_H
#define VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/Proxy/EpilogFwd.h"

#include "particle_variables.h"

/**
 * @namespace vars
 * @brief Namespace for organizing generic variables which act on interactions.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions. Each variable is implemented as a function which takes
 * an interaction object as an argument and returns a double. The function
 * should be templated on the type of interaction object if the variable is
 * intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the
 * pvars namespace, which is used for organizing variables which act on single
 * particles.
 */
namespace vars
{
    /**
     * @brief Variable for the true neutrino energy.
     * @units GeV
     * @tparam T the type of object.
     * @param obj the interaction to apply the variable on.
     * @return the true neutrino energy.
    */
    template<class T>
        double true_neutrino_energy(const T & obj) { return obj.nu_energy_init; }

    /**
     * @brief Variable for the true neutrino baseline. Currently, if the
     * distance is not available, the value is set to 585.0 meters.
     * @units meters
     * @tparam T the type of object.
     * @param obj interaction to apply the variable on.
     * @return the true neutrino baseline.
    */
    template<class T>
        double true_neutrino_baseline(const T & obj) { return obj.nu_distance_travel > 0 ? double(obj.nu_distance_travel) : 585.0; }

    /**
     * @brief Variable for the true neutrino PDG code.
     * @units none
     * @tparam T the type of object.
     * @param obj interaction to apply the variable on.
     * @return the true neutrino PDG code.
    */
    template<class T>
        double true_neutrino_pdg(const T & obj) { return obj.nu_pdg_code; }

    /**
     * @brief Variable for the true neutrino current value.
     * @units none
     * @tparam T the type of object.
     * @param obj interaction to apply the variable on.
     * @return the true neutrino current value.
    */
    template<class T>
        double true_neutrino_cc(const T & obj) { return obj.nu_current_type; }

    /**
     * @brief Variable for total visible energy of interaction.
     * @units GeV
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
    */
    template<class T>
        double visible_energy(const T & interaction)
        {
            double energy(0);
            for(const auto & p : interaction.particles)
            {
                if(p.is_primary)
                {
                    energy = particle::energy(p);
                    if(p.pid == 2) energy += MUON_MASS;
                    else if(p.pid == 3) energy += PION_MASS;
                }
            }
	        return energy/1000.0;
        }
}
#endif // VARIABLES_H
