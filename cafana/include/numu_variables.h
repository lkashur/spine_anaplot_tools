/**
 * @file numu_variables.h
 * @brief Header file for definitions of selection variables in the context
 * of the numu analyses.
 * @author mueller@fnal.gov
*/
#ifndef NUMU_VARIABLES_H
#define NUMU_VARIABLES_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "variables.h"
#include "cuts.h"

#include <algorithm>
#include <cmath>

namespace vars
{
   /**
     * @brief Finds the index corresponding to the leading particle of the specifed
     * particle type.
     * @tparam T the type of intearction (true or reco).
     * @param interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
    */
    template <class T>
        size_t leading_particle_index(const T & interaction, uint16_t pid)
        {
            double leading_ke(0);
            size_t index(0);
            for(size_t i(0); i < interaction.particles.size(); ++i)
            {
                const auto & p = interaction.particles[i];
                double energy(p.csda_ke);
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = ke_init(p);
                if(p.pid == pid && energy > leading_ke)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            return index;
        }

    /**
     * @brief Variable for the softmax score of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the softmax score of the leading muon.
    */
    template<class T>
        double leading_muon_softmax(const T & interaction)
        {
            auto & m(interaction.particles[leading_particle_index(interaction, 2)]);
            return m.pid_scores[2];
        }
    
    /**
     * @brief Variable for the softmax score of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the softmax score of the leading proton.
    */
    template<class T>
        double leading_proton_softmax(const T & interaction)
        {
            auto & p(interaction.particles[leading_particle_index(interaction, 4)]);
            return p.pid_scores[4];
        }

    /**
     * Variable for the transverse momentum of a particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the transverse momentum of the particle
    */
    template<class T>
        double transverse_momentum(const T & particle)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::sqrt(std::pow(particle.truth_momentum[0], 2) + std::pow(particle.truth_momentum[1], 2));
            else
                return std::sqrt(std::pow(particle.momentum[0], 2) + std::pow(particle.momentum[1], 2));
        }

    /**
     * Variable for the polar angle (w.r.t the z-axis) of the particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the polar angle of the particle.
    */
    template<class T>
        double polar_angle(const T & particle)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(particle.truth_start_dir[2]);
            else
                return std::acos(particle.start_dir[2]);
        }

    /**
     * Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the azimuthal angle of the particle.
    */
    template<class T>
        double azimuthal_angle(const T & particle)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(particle.truth_start_dir[0] / std::sqrt(std::pow(particle.truth_start_dir[0], 2) + std::pow(particle.truth_start_dir[1], 2)));
            else
                return std::acos(particle.start_dir[0] / std::sqrt(std::pow(particle.start_dir[0], 2) + std::pow(particle.start_dir[1], 2)));
        }

    /**
     * Variable for finding the leading muon kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    template<class T>
        double leading_muon_ke(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            double energy(csda_ke(interaction.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = ke_init(interaction.particles[i]);
            return energy;
        }
    
    /**
     * Variable for finding the leading proton kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    template<class T>
        double leading_proton_ke(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            double energy(csda_ke(interaction.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = ke_init(interaction.particles[i]);
            return energy;
        }

    /**
     * Variable for the transverse momentum of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading muon.
    */
    template<class T>
        double leading_muon_pt(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return transverse_momentum(interaction.particles[i]);
        }

    /**
     * Variable for the transverse momentum of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading proton.
    */
    template<class T>
        double leading_proton_pt(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return transverse_momentum(interaction.particles[i]);
        }

    /**
     * Variable for the muon polar angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the polar angle of the leading muon.
    */
    template<class T>
        double muon_polar_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return polar_angle(interaction.particles[i]);
        }

    /**
     * Variable for the muon azimuthal angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the azimuthal angle of the leading muon.
     */
    template<class T>
        double muon_azimuthal_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return azimuthal_angle(interaction.particles[i]);
        }

    /**
     * Variable for the opening angle between leading muon and
     * proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
    */
    template<class T>
        double opening_angle(const T & interaction)
        {
            auto & m(interaction.particles[leading_particle_index(interaction, 2)]);
            auto & p(interaction.particles[leading_particle_index(interaction, 4)]);
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return std::acos(m.truth_start_dir[0] * p.truth_start_dir[0] + m.truth_start_dir[1] * p.truth_start_dir[1] + m.truth_start_dir[2] * p.truth_start_dir[2]);
            else
                return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
    
    /**
     * Variable for the transverse momentum of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
    */
    template<class T>
        double interaction_pt(const T & interaction)
        {
            double px(0), py(0);
            for(const auto & p : interaction.particles)
                if(p.is_primary)
                {
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        px += p.truth_momentum[0];
                        py += p.truth_momentum[1];
                    }
                    else
                    {
                        px += p.momentum[0];
                        py += p.momentum[1];
                    }
                }
            return std::sqrt(std::pow(px, 2) + std::pow(py, 2));
        }

    /**
     * Variable for phi_T of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the phi_T of the interaction.
    */
    template<class T>
        double phiT(const T & interaction)
        {
            double lpx(0), lpy(0), hpx(0), hpy(0);
            for(const auto & p : interaction.particles)
                if(cuts::final_state_signal(p))
                {
                    if(p.pid > 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            hpx += p.truth_momentum[0];
                            hpy += p.truth_momentum[1];
                        }
                        else
                        {
                            hpx += p.momentum[0];
                            hpy += p.momentum[1];
                        }
                    }
                    else if(p.pid == 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            lpx += p.truth_momentum[0];
                            lpy += p.truth_momentum[1];
                        }
                        else
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                    }
                }
            return std::acos((-hpx * lpx - hpy * lpy) / (std::sqrt(std::pow(hpx, 2) + std::pow(hpy, 2)) * std::sqrt(std::pow(lpx, 2) + std::pow(lpy, 2))));
        }

    /**
     * Variable for alpha_T of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the phi_T of the interaction.
    */
    template<class T>
        double alphaT(const T & interaction)
        {
            double lpx(0), lpy(0), px(0), py(0);
            for(const auto & p : interaction.particles)
                if(cuts::final_state_signal(p))
                {
                    if(p.pid <= 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            lpx += p.truth_momentum[0];
                            lpy += p.truth_momentum[1];
                        }
                        else
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                    }
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        px += p.truth_momentum[0];
                        py += p.truth_momentum[1];
                    }
                    else
                    {
                        px += p.momentum[0];
                        py += p.momentum[1];
                    }
                }
            return std::acos((-px * lpx - py * lpy) / (std::sqrt(std::pow(px, 2) + std::pow(py, 2)) * std::sqrt(std::pow(lpx, 2) + std::pow(lpy, 2))));
        }

    /**
     * Variable for the muon softmax score for the leading muon of the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the muon softmax score of the leading muon.
    */
    template<class T>
        double muon_softmax(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return interaction.particles[i].pid_scores[2];
        }

    /**
     * Variable for the proton softmax score for the leading proton of the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the proton softmax score of the leading proton.
    */
    template<class T>
        double proton_softmax(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return interaction.particles[i].pid_scores[4];
        }
}
#endif