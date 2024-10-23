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
}
#endif // PARTICLE_VARIABLES_H