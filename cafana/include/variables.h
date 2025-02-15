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

#include "include/particle_variables.h"
#include "include/particle_cuts.h"
#include "include/cuts.h"
#include "include/utilities.h"
#include "include/particle_utilities.h"

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
     * @brief Variable for the neutrino ID of the interaction.
     * @details This variable is intended to provide a unique identifier for
     * each parent neutrino within the event record. This number is assigned
     * starting at 0 for the first neutrino in the event and is incremented
     * for each subsequent neutrino. Non-neutrino interactions are assigned
     * a value of -1.
     * @param obj the interaction to apply the variable on.
     * @return the neutrino ID.
     */
    double neutrino_id(const caf::SRInteractionTruthDLPProxy & obj) { return obj.nu_id; }

    /**
     * @brief Variable for the containment status of the interaction.
     * @details The containment status is determined upstream in the SPINE
     * reconstruction and is based on the set of all points in the interaction,
     * which must be contained within the volume of the TPC that created them.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the containment status of the interaction.
     */
    template<class T>
        double containment(const T & obj) { return obj.is_contained; }

    /**
     * @brief Variable for the fiducial volume status of the interaction.
     * @details The fiducial volume status is determined upstream in the SPINE
     * reconstruction and is a requirement that the interaction vertex is within
     * the fiducial volume of the TPC.
     */
    template<class T>
        double fiducial(const T & obj) { return obj.is_fiducial; }

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
                if(pcuts::final_state_signal(p))
                {
                    energy += pvars::energy(p);
                    if(p.pid == 4) energy -= PROTON_MASS - PROTON_BINDING_ENERGY;
                }
            }
            return energy/1000.0;
        }

    /**
     * @brief Variable for total visible energy of interaction, including
     * sub-threshold particles.
     * @details This function calculates the total visible energy of the
     * interaction by summing the energy of all particles that are identified
     * as counting towards the final state of the interaction. Sub-threshold
     * particles are included calorimetrically.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
        double visible_energy_calosub(const T & obj)
        {
            double energy(0);
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    energy += pvars::energy(p);
                    if(p.pid == 4) energy -= PROTON_MASS - PROTON_BINDING_ENERGY;
                }
                else if(pcuts::is_primary(p))
                    energy += p.calo_ke;
            }
            return energy/1000.0;
        }

    /**
     * @brief Variable for the flash time of the interaction.
     * @details The flash time is the time of the flash observed in the PMTs
     * and associated with the charge deposition in the interaction using the
     * OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash time of the interaction.
     */
    template<class T>
        double flash_time(const T & obj)
        {
            if(obj.flash_times.size() > 0)
                return obj.flash_times[0];
            return PLACEHOLDERVALUE;
        }

    /**
     * @brief Variable for the flash total photoelectron count of the
     * interaction.
     * @details The flash total photoelectron count is the total number of
     * photoelectrons observed in the PMTs and associated with the charge
     * deposition in the interaction using the OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash total photoelectron count of the interaction.
     */
    template<class T>
        double flash_total_pe(const T & obj) { return obj.flash_total_pe; }

    /**
     * @brief Variable for the flash hypothesis total photoelectron count of
     * the interaction.
     * @details The flash hypothesis total photoelectron count is the total
     * number of photoelectrons predicted by OpT0Finder for the interaction in
     * the flash associated with the charge deposition in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash hypothesis total photoelectron count of the interaction.
     */
    template<class T>
        double flash_hypothesis(const T & obj) { return obj.flash_hypo_pe; }

    /**
     * @brief Variable for the truth x-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @param obj the interaction to apply the variable on.
     * @return the truth x-coordinate of the interaction vertex.
     */
    double truth_vertex_x(const caf::SRInteractionTruthDLPProxy & obj) { return obj.vertex[0]; }

    /**
     * @brief Variable for the truth y-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @param obj the interaction to apply the variable on.
     * @return the truth y-coordinate of the interaction vertex.
     */
    double truth_vertex_y(const caf::SRInteractionTruthDLPProxy & obj) { return obj.vertex[1]; }

    /**
     * @brief Variable for the truth z-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @param obj the interaction to apply the variable on.
     * @return the truth z-coordinate of the interaction vertex.
     */
    double truth_vertex_z(const caf::SRInteractionTruthDLPProxy & obj) { return obj.vertex[2]; }

    /**
     * @brief Variable for the x-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the x-coordinate of the interaction vertex.
     */
    template<class T>
        double vertex_x(const T & obj) { return obj.vertex[0]; }

    /**
     * @brief Variable for the y-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the y-coordinate of the interaction vertex.
     */
    template<class T>
        double vertex_y(const T & obj) { return obj.vertex[1]; }

    /**
     * @brief Variable for the z-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the z-coordinate of the interaction vertex.
     */
    template<class T>
        double vertex_z(const T & obj) { return obj.vertex[2]; }

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only particles identified as contributing to the final state.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of all particles that are
     * identified as counting towards the final state of the interaction. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double dpT(const T & obj)
        {
            utilities::three_vector pt = {0, 0, 0};
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    // Sum up the transverse momentum of all final state particles
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                    pt = utilities::add(pt, this_pt);
                }
            }
            return utilities::magnitude(pt);
        }

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only the leading charged lepton and proton.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of the leading charged
     * lepton and proton. The neutrino direction is assumed to either be the
     * BNB axis direction (z-axis) or the unit vector pointing from the NuMI
     * target to the interaction vertex. See @ref utilities::transverse_momentum
     * for details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the leading charged lepton and proton.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double dpT_lp(const T & obj)
        {
            
            utilities::three_vector l_pt = {0, 0, 0};
            utilities::three_vector p_pt = {0, 0, 0};
            double l_ke(0), p_ke(0);
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    // Find the leading charged lepton and proton
                    if((p.pid == 1 || p.pid == 2) && pvars::ke(p) > l_ke)
                    {
                        l_ke = pvars::ke(p);
                        utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                        utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                        l_pt = utilities::transverse_momentum(momentum, vtx);
                    }
                    else if(p.pid == 4 && pvars::ke(p) > p_ke)
                    {
                        p_ke = pvars::ke(p);
                        utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                        utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                        p_pt = utilities::transverse_momentum(momentum, vtx);
                    }
                }
            }
            if(l_ke == 0 || p_ke == 0)
                return PLACEHOLDERVALUE;
            else
                return utilities::magnitude(utilities::add(l_pt, p_pt));
        }

    /**
     * @brief Variable for phi_T of the interaction.
     * @details phi_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the leading muon and the total hadronic
     * system. This variable is sensitive to the presence of F.S.I. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the phi_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double phiT(const T & obj)
        {
            utilities::three_vector lepton_pt = {0, 0, 0};
            utilities::three_vector hadronic_pt = {0, 0, 0};
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    // There should only be one lepton, so replace the lepton
                    // transverse momentum if the particle is a lepton.
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                    if(p.pid == 1 || p.pid == 2)
                        lepton_pt = this_pt;
                    // The total hadronic system is treated as a single object.
                    else if(p.pid > 2)
                        hadronic_pt = utilities::add(hadronic_pt, this_pt);
                }
            }
            return std::acos(-1 * utilities::dot_product(lepton_pt, hadronic_pt) / (utilities::magnitude(lepton_pt) * utilities::magnitude(hadronic_pt)));
        }

    /**
     * @brief Variable for alpha_T of the interaction.
     * @details alpha_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the total hadronic system and the
     * outgoing lepton. The neutrino direction is assumed to either be the BNB
     * axis direction (z-axis) or the unit vector pointing from the NuMI target
     * to the interaction vertex. See @ref utilities::transverse_momentum for
     * details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the alpha_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double alphaT(const T & obj)
        {
            utilities::three_vector lepton_pt = {0, 0, 0};
            utilities::three_vector total_pt = {0, 0, 0};
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    // There should only be one lepton, so replace the lepton
                    // transverse momentum if the particle is a lepton.
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                    if(p.pid == 1 || p.pid == 2)
                        lepton_pt = this_pt;
                    total_pt = utilities::add(total_pt, this_pt);
                }
            }
            return std::acos(-1 * utilities::dot_product(total_pt, lepton_pt) / (utilities::magnitude(total_pt) * utilities::magnitude(lepton_pt)));
        }

    /**
     * @brief Variable for the missing longitudinal momentum of the
     * interaction.
     * @details The missing longitudinal momentum is calculated as the
     * difference between the total longitudinal momentum of the final state
     * particles and the best estimate of the neutrino energy. The neutrino
     * energy is calculated using @ref vars::visible_energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the missing longitudinal momentum of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double dpL(const T & obj)
        {
            utilities::three_vector lepton_pl = {0, 0, 0};
            utilities::three_vector hadronic_pl = {0, 0, 0};
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    // There should only be one lepton, so replace the lepton
                    // transverse momentum if the particle is a lepton.
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    utilities::three_vector this_pl = utilities::longitudinal_momentum(momentum, vtx);
                    if(p.pid == 1 || p.pid == 2)
                        lepton_pl = this_pl;
                    // The total hadronic system is treated as a single object.
                    else if(p.pid > 2)
                        hadronic_pl = utilities::add(hadronic_pl, this_pl);
                }
            }
            return utilities::magnitude(utilities::add(hadronic_pl, lepton_pl)) - 1000*vars::visible_energy(obj);
        }

    /**
     * @brief Variable for the missing longitudinal momentum of the interaction
     * counting only the leading charged lepton and proton.
     * @details The missing longitudinal momentum is calculated as the
     * difference between the total longitudinal momentum of the leading charged
     * lepton and proton and the best estimate of the neutrino energy. The
     * neutrino energy is calculated using @ref vars::visible_energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the missing longitudinal momentum of the leading charged lepton
     * and proton.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double dpL_lp(const T & obj)
        {
            utilities::three_vector l_pl = {0, 0, 0};
            utilities::three_vector p_pl = {0, 0, 0};
            double l_ke(0), p_ke(0);
            for(const auto & p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                {
                    // Find the leading charged lepton and proton
                    if((p.pid == 1 || p.pid == 2) && pvars::ke(p) > l_ke)
                    {
                        l_ke = pvars::ke(p);
                        utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                        utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                        l_pl = utilities::longitudinal_momentum(momentum, vtx);
                    }
                    else if(p.pid == 4 && pvars::ke(p) > p_ke)
                    {
                        p_ke = pvars::ke(p);
                        utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                        utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                        p_pl = utilities::longitudinal_momentum(momentum, vtx);
                    }
                }
            }
            if(l_ke == 0 || p_ke == 0)
                return PLACEHOLDERVALUE;
            else
                return utilities::magnitude(utilities::add(l_pl, p_pl)) - 1000*vars::visible_energy(obj);
        }

    /**
     * @brief Variable for the estimate of the momentum of the struck nucleon.
     * @details The estimate of the momentum of the struck nucleon is calculated
     * as the quadrature sum of the transverse momentum (see @ref vars::dpT) and
     * the missing longitudinal momentum (see @ref vars::dpL).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the estimate of the momentum of the struck nucleon.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double pn(const T & obj) { return std::sqrt(std::pow(vars::dpT(obj), 2) + std::pow(vars::dpL(obj), 2)); }

    /**
     * @brief Variable for the estimate of the momentum of the struck nucleon
     * counting only the leading charged lepton and proton.
     * @details The estimate of the momentum of the struck nucleon is calculated
     * as the quadrature sum of the transverse momentum (see @ref vars::dpT_lp)
     * and the missing longitudinal momentum (see @ref vars::dpL_lp).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the estimate of the momentum of the struck nucleon.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double pn_lp(const T & obj) { return std::sqrt(std::pow(vars::dpT_lp(obj), 2) + std::pow(vars::dpL_lp(obj), 2)); }
}
#endif // VARIABLES_H
