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
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"

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
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 0: 1mu1p (contained and fiducial)
     * 1: 1mu1p (not contained or not fiducial)
     * 2: 1muNp (N > 1, contained and fiducial)
     * 3: 1muNp (N > 1, not contained or fiducial)
     * 4: 1muX (not 1muNp, contained and fiducial)
     * 5: 1muX (not 1muNp, not contained or fiducial)
     * 6: Other CC nu
     * 7: Other NC nu
     * 8: Cosmic
     * @param obj The interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    double category(const caf::SRInteractionTruthDLPProxy & obj)
    {
        double cat(8);
        if(cuts::muon2024::signal_1mu1p(obj)) cat = 0;
        else if(cuts::muon2024::nonsignal_1mu1p(obj)) cat = 1;
        else if(cuts::muon2024::signal_1muNp(obj)) cat = 2;
        else if(cuts::muon2024::nonsignal_1muNp(obj)) cat = 3;
        else if(cuts::muon2024::signal_1muX(obj)) cat = 4;
        else if(cuts::muon2024::nonsignal_1muX(obj)) cat = 5;
        else if(cuts::neutrino(obj) && cuts::iscc(obj)) cat = 6;
        else if(cuts::neutrino(obj) && !cuts::iscc(obj)) cat = 7;
        return cat;
    }

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
            return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
}
#endif // VARS_MUON2024_H