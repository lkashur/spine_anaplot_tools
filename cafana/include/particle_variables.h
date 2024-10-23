/**
 * @file particle_variables.h
 * @brief Header file for definitions of variables which act on single
 * particles.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_VARIABLES_H
#define PARTICLE_VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

namespace vars::particle
{
    /**
     * @brief Variable for the best estimate of the particle energy. At the
     * most basic decision level, this is based on the shower/track
     * designation. Showers can only be reconstructed calorimetrically, while
     * tracks can be reconstructed calorimetrically, by range (if contained),
     * or by multiple scattering (if exiting).
     * @tparam T the type of particle.
     * @param p the particle to apply the variable on.
     * @return the best estimate of the particle energy.
     */
    template<class T>
        double energy(const T & p)
        {
            double energy = 0;
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
            {
                energy += p.energy_deposit;
            }
            else
            {
                // Check if the particle is a shower.
                if(p.pid < 2) energy += p.calo_ke;
                else
                {
                    if(p.is_contained) energy += p.csda_ke;
                    else energy += p.mcs_ke;
                }
            }
            return energy;
        }
    
    /**
     * @brief Variable for the transverse momentum of a particle.
     * @details The transverse momentum is defined as the square root of the
     * sum of the squares of the x and y components of the momentum. This
     * variable is useful for identifying particles which are produced in a
     * transverse direction to the beam.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the transverse momentum of the particle.
     */
    template<class T>
        double transverse_momentum(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::sqrt(std::pow(p.truth_momentum[0], 2) + std::pow(p.truth_momentum[1], 2));
            else
                return std::sqrt(std::pow(p.momentum[0], 2) + std::pow(p.momentum[1], 2));
        }

    /**
     * @brief Variable for the polar angle (w.r.t the z-axis) of the particle.
     * @details The polar angle is defined as the arccosine of the z-component
     * of the momentum vector. This variable is useful for identifying particles
     * which are produced transversely to the beam.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the polar angle of the particle.
     */
    template<class T>
        double polar_angle(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(p.truth_start_dir[2]);
            else
                return std::acos(p.start_dir[2]);
        }

    /**
     * @brief Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @details The azimuthal angle is defined as the arccosine of the x-component
     * of the momentum vector divided by the square root of the sum of the squares
     * of the x and y components of the momentum vector.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the azimuthal angle of the particle.
     */
    template<class T>
        double azimuthal_angle(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(p.truth_start_dir[0] / std::sqrt(std::pow(p.truth_start_dir[0], 2) + std::pow(p.truth_start_dir[1], 2)));
            else
                return std::acos(p.start_dir[0] / std::sqrt(std::pow(p.start_dir[0], 2) + std::pow(p.start_dir[1], 2)));
        }
}
#endif // PARTICLE_VARIABLES_H