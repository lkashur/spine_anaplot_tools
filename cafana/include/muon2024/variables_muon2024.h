/**
 * @file vars_muon2024.h
 * @brief Header file for definitions of analysis variables specific to the
 * muon2024 analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions specific to the muon2024
 * analysis. Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author mueller@fnal.gov
 */
#ifndef VARS_MUON2024_H
#define VARS_MUON2024_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/utilities.h"

/**
 * @namespace vars::muon2024
 * @brief Namespace for organizing variables specific to the muon2024 analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the muon2024 analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::muon2024
{
    /**
     * @brief Variable for the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
     */
    template<class T>
        double opening_angle(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return std::acos(m.truth_start_dir[0] * p.truth_start_dir[0] + m.truth_start_dir[1] * p.truth_start_dir[1] + m.truth_start_dir[2] * p.truth_start_dir[2]);
            else
                return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
}
#endif // VARS_MUON2024_H