/**
 * @file particle_cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select particles. Each cut is implemented as a function which takes an
 * particle object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections/interaction cuts.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_CUTS_H
#define PARTICLE_CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "particle_variables.h"

/**
 * @namespace pcuts
 * @brief Namespace for organizing generic cuts which act on particles.
 * @details This namespace is intended to be used for organizing cuts which act
 * on particles. Each cut is implemented as a function which takes a
 * particle object as an argument and returns a boolean. The function should
 * be templated on the type of particle object if the cut is intended to be
 * used on both true and reconstructed particles.
 */
namespace pcuts
{
    /**
     * @brief Check if the particle is a primary particle.
     * @details This function checks if the particle is a primary particle.
     * Primary designation is handled upstream in SPINE and is based on the
     * softmax primary/secondary scores that are assigned to each particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a primary particle.
     */
    template<class T>
        bool is_primary(const T & p)
        {
            return p.is_primary;
        }
    
    /**
     * @brief Check if particle is contained.
     * @details This function checks if the particle is contained.
     * Containment designation is handled upstream in SPINE and considers
     * a 5 cm buffer with regard to the boundaries of the detector.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is contained.
     *
     */
    template<class T>
        bool is_contained(const T & p)
        {
	    return p.is_contained;
        }

    /**
     * @brief Check if the particle meets final state signal requirements.
     * @details must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(is_primary(p))
            {
                double energy(pvars::ke(p));
                if((PIDFUNC(p) == 2 && energy > 143.425) || (PIDFUNC(p) != 2 && PIDFUNC(p) < 4 && energy > 25) || (PIDFUNC(p) == 4 && energy > 50))
                    passes = true;
            }
            return passes;
        }

    /**
     * @brief Check if the particle is throughgoing.
     * @details This function checks if the particle is throughgoing. A
     * throughgoing particle is defined as a particle which has both ends
     * of the track near the boundary of the detector. This is only applicable
     * to tracks as it is somewhat nonsensical for showers.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is throughgoing.
     */
    template<class T>
        bool throughgoing(const T & p)
        {
            return PIDFUNC(p) > 1 && utilities::near_boundary(p.start) && utilities::near_boundary(p.end);
        }
}
#endif // PARTICLE_CUTS_H
