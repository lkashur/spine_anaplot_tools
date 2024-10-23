/**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_H
#define CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "particle_variables.h"

/**
 * @namespace cuts
 * @brief Namespace for organizing generic cuts which act on interactions.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. The function should
 * be templated on the type of interaction object if the cut is intended to be
 * used on both true and reconstructed interactions.
 */
namespace cuts
{
    /**
     * @brief Check if the particle meets final state signal requirements.
     * Particles must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param particle to check.
     * @return true if the particle is a final state signal particle.
    */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(p.is_primary)
            {
                double energy(vars::particle::energy(p));
                if((p.pid == 2 && energy > 143.425) || (p.pid != 2 && p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 50))
                    passes = true;
            }
            return passes;
        }

    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
    */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & interaction)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : interaction.particles)
            {
                if(final_state_signal(p))
                    ++counts[p.pid];
            }
            return counts;
        }
    
    /**
     * @brief Apply a cut on the validity of the flash match.
     * @tparam T the type of interaction (true or reco).
     * @param interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
    */
    template<class T>
        bool valid_flashmatch(const T & interaction) { return !std::isnan(interaction.flash_time) && interaction.fmatched == 1; }

    /**
     * @brief Apply no cut (all interactions/particles passed).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to select on.
     * @return true (always).
    */
    template<class T>
        bool no_cut(const T & obj) { return true; }

    /**
     * @brief Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
    */
    template<class T>
        bool fiducial_cut(const T & interaction)
        {
            return interaction.is_fiducial && !(interaction.vertex[0] > 210.215 && interaction.vertex[1] > 60 && (interaction.vertex[2] > 290 && interaction.vertex[2] < 390));
        }
    
    /**
     * @brief Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is contained.
    */
    template<class T>
        bool containment_cut(const T & interaction) { return interaction.is_contained; }

    /**
     * @brief Apply a 1mu1p topological cut. The interaction must have a topology
     * matching 1mu1p as defined by the conditions in the topology() function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1mu1p topology.
    */
    template<class T>
        bool topological_1mu1p_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] == 1;
        }
    
    /**
     * @brief Apply a 1muNp topological cut. The interaction must have a topology
     * matching 1muNp as defined by the conditions in the count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muNp topology.
    */
    template<class T>
        bool topological_1muNp_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] >= 1;
        }

    /**
     * @brief Apply a 1muX topological cut. The interaction must have a topology
     * matching 1muX as defined by the conditions in the count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muX topology.
    */
    template<class T>
        bool topological_1muX_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[2] == 1;
        }
    
    /**
     * @brief Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
    */
    template<class T>
        bool flash_cut(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6);
        }

    /**
     * @brief Apply a fiducial and containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
    */
    template<class T>
        bool fiducial_containment_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * @brief Apply a fiducial, containment, and topological (1mu1p) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
    */
    template<class T>
        bool fiducial_containment_topological_1mu1p_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1mu1p_cut<T>(interaction); }

    /**
     * @brief Apply a fiducial, containment, and topological (1muNp) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
    */
    template<class T>
        bool fiducial_containment_topological_1muNp_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1muNp_cut<T>(interaction); }

    /**
     * @brief Apply a fiducial, containment, and topological (1muX) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
    */
    template<class T>
        bool fiducial_containment_topological_1muX_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1muX_cut<T>(interaction); }

    /**
     * @brief Apply a fiducial, containment, topological (1mu1p), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
    */
    template<class T>
        bool all_1mu1p_cut(const T & interaction) { return topological_1mu1p_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * @brief Apply a fiducial, containment, topological (1muNp), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
    */
    template<class T>
        bool all_1muNp_cut(const T & interaction) { return topological_1muNp_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }
    
    /**
     * @brief Apply a fiducial, containment, topological (1muX), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
    */
}
#endif
