/**
 * @file particle_variables.h
 * @brief Header file for definitions of variables which act on single
 * particles.
 * @details This file contains definitions of variables which act on single
 * particles. Each variable is implemented as a function which takes a particle
 * object as an argument and returns a double. These variables are intended to
 * be used to define more complex variables which act on interactions.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_VARIABLES_H
#define PARTICLE_VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include "include/particle_utilities.h"

/**
 * @namespace pvars
 * @brief Namespace for organizing generic variables which act on single
 * particles.
 * @details This namespace is intended to be used for organizing variables which
 * act on single particles. Each variable is implemented as a function which
 * takes a particle object as an argument and returns a double. The function
 * should be templated on the type of particle object if the variable is
 * intended to be used on both true and reconstructed particles.
 * @note The namespace is intended to be used in conjunction with the
 * vars namespace, which is used for organizing variables which act on
 * interactions.
 */
namespace pvars
{
    /**
     * @brief Variable for the particle's primary classification.
     * @details This variable returns the primary classification of the particle.
     * The primary classification is determined upstream in the SPINE
     * reconstruction and is based on the softmax scores of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the primary classification of the particle.
     */
    template<class T>
        double primary_classification(const T & p)
        {
            return p.is_primary ? 1 : 0;
        }

    /**
     * @brief Variable for assigning primary classification based on the
     * particle's softmax scores.
     * @details This variable assigns a primary classification based on the
     * softmax scores of the particle. This function places a relaxed threshold
     * on the primary softmax score to reduce observed inefficiencies in the
     * primary classification.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the primary classification of the particle.
     */
    template<class T>
        double custom_primary_classification(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            {
                return p.is_primary ? 1 : 0;
            }
            else
            {
                return p.primary_scores[1] > 0.10 ? 1 : 0;
            }
        }

    /**
     * @brief Variable for the particle's PID.
     * @details This variable returns the PID of the particle. The PID is
     * determined by the softmax scores of the particle. This function uses the
     * "nominal" PID decision that is made upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
        double pid(const T & p)
        {
            return p.pid;
        }

    /**
     * @brief Variable for assigning PID based on the particle's softmax scores.
     * @details This variable assigns a PID based on the softmax scores of the
     * particle. Nominally, the PID is assigned based on the highest softmax
     * score, but the PID can be overridden directly by this function.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
        double custom_pid(const T & p)
        {
            double pid = std::numeric_limits<double>::quiet_NaN();
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            {
                pid = p.pid;
            }
            else
            {
                if(p.pid_scores[2] > 0.25)
                    pid = 2;
                else
                {
                    size_t high_index(0);
                    for(size_t i(0); i < 5; ++i)
                        if(p.pid_scores[i] > p.pid_scores[high_index]) high_index = i;
                    pid = high_index;
                }
            }
            return pid;
        }

    /**
     * @brief Variable for the semantic type of the particle.
     * @details This variable returns the semantic type of the particle. The
     * semantic type is determined by majority-vote of the pixel-level semantic
     * segmentation of the particle. The semantic types are defined as follows:
     * 0: shower, 1: track, 2: Michel electron, 3: delta electron,
     * 4: low-energy, 5: ghost, and -1: unknown.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the semantic type of the particle.
     */
    template<class T>
        double semantic_type(const T & p)
        {
            return p.shape;
        }

    /**
     * @brief Variable for the best-match IoU of the particle.
     * @details The best-match IoU is the intersection over union of the
     * points belonging to a pair of reconstructed and true particles. The
     * best-match IoU is calculated upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the best-match IoU of the particle.
     */
    template<class T>
        double iou(const T & p)
        {
            if(p.match_ids.size() > 0)
                return p.match_overlaps[0];
            else 
                return PLACEHOLDERVALUE;
        }

    /**
     * @brief Variable for the mass of the particle.
     * @details The mass of the particle is determined by the PID of the
     * particle. This couples the PID to the mass of the particle, so it is
     * necessary to use the appropriate PID function rather than the in-built
     * PID attribute.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the mass of the particle.
     */
    template<class T>
        double mass(const T & p)
        {
            double mass(0);
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            {
                mass = p.mass;
            }
            else
            {
                switch(int(PIDFUNC(p)))
                {
                    case 0:
                        mass = 0;
                        break;
                    case 1:
                        mass = ELECTRON_MASS;
                        break;
                    case 2:
                        mass = MUON_MASS;
                        break;
                    case 3:
                        mass = PION_MASS;
                        break;
                    case 4:
                        mass = PROTON_MASS;
                        break;
                    default:
                        mass = PLACEHOLDERVALUE;
                        break;
                }
            }
            return mass;
        }

    /**
     * @brief Variable for the CSDA kinetic energy of the particle.
     * @details The CSDA kinetic energy is calculated upstream in the SPINE
     * reconstruction by relating the length of the track to the energy loss
     * of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the CSDA kinetic energy of the particle.
     * @note This function is only valid for particles which are contained in
     * the detector.
     */
    template<class T>
        double csda_ke(const T & p)
        {
            return p.csda_ke_per_pid[PIDFUNC(p)];
        }

    /**
     * @brief Variable for the MCS kinetic energy of the particle.
     * @details The MCS kinetic energy is calculated upstream in the SPINE
     * reconstruction by comparing the distribution of successive angles of
     * track segments to the expected distribution for multiple scattering.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the MCS kinetic energy of the particle.
     */
    template<class T>
        double mcs_ke(const T & p)
        {
            return p.mcs_ke_per_pid[PIDFUNC(p)];
        }

    /**
     * @brief Variable for the calorimetric kinetic energy of the particle.
     * @details The calorimetric kinetic energy is calculated upstream in the
     * SPINE reconstruction as the sum of the energy of each spacepoint in the
     * particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the calorimetric kinetic energy of the particle.
     */
    template<class T>
        double calo_ke(const T & p)
        {
            return p.calo_ke;
        }

    /**
     * @brief Variable for true particle starting kinetic energy.
     * @details The starting kinetic energy is defined as the total energy
     * minus the rest mass energy of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the starting kinetic energy of the particle.
     */
    template<class T>
        double ke(const T & p)
        {
            double energy(0);
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            {
                energy = p.energy_init - mass(p);
            }
            else
            {
                if(PIDFUNC(p) < 2) [[likely]]
                    energy += calo_ke(p);
                else
                {
                    if(p.is_contained) energy += csda_ke(p);
                    else energy += mcs_ke(p);
                }
            }
            return energy;
        }

    /**
     * @brief Variable for particle CSDA kinetic energy.
     * @details The kinetic energy as calculated upstream in SPINE
     * using the continuous slowing down approximation.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle CSDA kinetic energy.
     */
    /*
    template<class T>
        double csda_ke(const T & p)
        {
	    return p.csda_ke;
        }
    */

    /**
     * @brief Variable for particle calorimetric kinetic energy.
     * @details The kinetic energy as calculated upstream in SPINE
     * using the sum of calorimetric energy deposits.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle calorimetric kinetic energy.
     */
    /*
    template<class T>
        double calo_ke(const T & p)
        {
	    return p.calo_ke;
	}
    */

    /**
     * @brief Variable for particle calorimetric kinetic energy, pre-shower correction.
     * @details The kinetic energy as calculated upstream in SPINE
     * using the sum of calorimetric energy deposits, without applying a shower correction factor.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle calorimetric kinetic energy, pre-shower correction.
     */
    template<class T>
        double calo_ke_pre_corr(const T & p)
        {
	    return (1/1.2359) * p.calo_ke;
        }

    template<class T>
        double custom_calo_ke(const T & p)
        {
	  return (1/52.356) * (0.77) * (47.55) * (1/0.762) * p.calo_ke;
        }

    /**
     * @brief Variable for particle MCS kinetic energy.
     * @details The kinetic energy as calculated upstream in SPINE
     * using multiple coulomb scattering.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle MCS kinetic energy.
     */
    /*
    template<class T>
        double mcs_ke(const T & p)
        {
	    return p.mcs_ke;
        }
    */

    /**
     * @brief Variable for the best estimate of the particle energy.
     * @details At the most basic decision level, this is based on the
     * shower/track designation. Showers can only be reconstructed
     * calorimetrically, while tracks can be reconstructed calorimetrically,
     * by range (if contained), or by multiple scattering (if exiting).
     * @tparam T the type of particle.
     * @param p the particle to apply the variable on.
     * @return the best estimate of the particle energy.
     */
    template<class T>
        double energy(const T & p)
        {
            double energy = ke(p) + mass(p);
            return energy;
        }

    /**
     * @brief Variable for the length of the particle track.
     * @details The length of the track is calculated upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the length of the particle track.
     */
    template<class T>
        double length(const T & p)
        {
            return p.length;
        }

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * a horizontal anode plane.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of a horizontal anode plane (e.g. the front induction plane of ICARUS).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
        double theta_xw_horizontal(const T & p)
        {
            // Wire orientation
            double wx(0), wy(1), wz(0);
            return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
        }

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle +60 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of +60 degrees from
     * the horizontal (e.g. the second induction / collection plane of ICARUS).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
        double theta_xw_p60(const T & p)
        {
            // Wire orientation
            double wx(0), wy(0.5 * std::sqrt(3)), wz(0.5);
            return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
        }
    
    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle -60 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of -60 degrees from
     * the horizontal (e.g. the third induction / collection plane of ICARUS).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
        double theta_xw_m60(const T & p)
        {
            // Wire orientation
            double wx(0), wy(-0.5 * std::sqrt(3)), wz(0.5);
            return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
        }
    
    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * a vertical anode plane.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of a vertical anode plane (e.g. the collection plane of SBND).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
        double theta_xw_vertical(const T & p)
        {
            // Wire orientation
            double wx(0), wy(1), wz(0);
            return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
        }

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle +30 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of +30 degrees from
     * the horizontal (e.g. the first induction plane of SBND).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
        double theta_xw_p30(const T & p)
        {
            // Wire orientation
            double wx(0), wy(0.5), wz(0.5 * std::sqrt(3));
            return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
        }

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle -30 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of -30 degrees from
     * the horizontal (e.g. the second induction plane of SBND).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
        double theta_xw_m30(const T & p)
        {
            // Wire orientation
            double wx(0), wy(-0.5), wz(0.5 * std::sqrt(3));
            return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
        }

    /**
     * @brief Variable for the x-coordinate of the particle starting point.
     * @details The starting point is the point at which the particle is created
     * and is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-coordinate of the particle starting point.
     */
    template<class T>
        double start_x(const T & p)
        {
            return p.start_point[0];
        }

    /**
     * @brief Variable for the y-coordinate of the particle starting point.
     * @details The starting point is the point at which the particle is created
     * and is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-coordinate of the particle starting point.
     */
    template<class T>
        double start_y(const T & p)
        {
            return p.start_point[1];
        }

    /**
     * @brief Variable for the z-coordinate of the particle starting point.
     * @details The starting point is the point at which the particle is created
     * and is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-coordinate of the particle starting point.
     */
    template<class T>
        double start_z(const T & p)
        {
            return p.start_point[2];
        }

    /**
     * @brief Variable for the x-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-coordinate of the particle end point.
     */
    template<class T>
        double end_x(const T & p)
        {
            return p.end_point[0];
        }

    /**
     * @brief Variable for the y-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-coordinate of the particle end point.
     */
    template<class T>
        double end_y(const T & p)
        {
            return p.end_point[1];
        }
    
    /**
     * @brief Variable for the z-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-coordinate of the particle end point.
     */
    template<class T>
        double end_z(const T & p)
        {
            return p.end_point[2];
        }

    /**
     * @brief Variable for the x-component of the particle momentum.
     * @details The momentum is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-component of the particle momentum.
     */
    template<class T>
        double px(const T & p)
        {
            return p.momentum[0];
        }
    
    /**
     * @brief Variable for the y-component of the particle momentum.
     * @details The momentum is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-component of the particle momentum.
     */
    template<class T>
        double py(const T & p)
        {
            return p.momentum[1];
        }

    /**
     * @brief Variable for the z-component of the particle momentum.
     * @details The momentum is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-component of the particle momentum.
     */
    template<class T>
        double pz(const T & p)
        {
            return p.momentum[2];
        }
    
    /**
     * @brief Variable for the transverse momentum of a particle.
     * @details This function calculates the transverse momentum of the
     * particle with respect to the assumed neutrino direction. The neutrino
     * direction is assumed to either be the BNB axis direction (z-axis) or the
     * unit vector pointing from the NuMI target to the interaction vertex.
     * See @ref utilities::transverse_momentum for details on the extraction of
     * the transverse momentum.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the transverse momentum of the particle.
     */
    template<class T>
        double dpT(const T & p)
        {
            utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
            utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
            utilities::three_vector pt = utilities::transverse_momentum(momentum, vtx);
            return utilities::magnitude(pt);
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
            return std::acos(p.start_dir[0] / std::sqrt(std::pow(p.start_dir[0], 2) + std::pow(p.start_dir[1], 2)));
        }

    /**
     * @brief Variable for the photon softmax score of the particle.
     * @details The photon softmax score represents the confidence that the
     * network has in the particle being a photon. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a photon.
     * @param p the particle to apply the variable on.
     * @return the photon softmax score of the particle.
     */
    double photon_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[0];
    }

    /**
     * @brief Variable for the electron softmax score of the particle.
     * @details The electron softmax score represents the confidence that the
     * network has in the particle being an electron. The score is between 0 and 1,
     * with 1 being the most confident that the particle is an electron.
     * @param p the particle to apply the variable on.
     * @return the electron softmax score of the particle.
     */
    double electron_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[1];
    }
    
    /**
     * @brief Variable for the muon softmax score of the particle.
     * @details The muon softmax score represents the confidence that the
     * network has in the particle being a muon. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a muon.
     * @param p the particle to apply the variable on.
     * @return the muon softmax score of the particle.
     */
    double muon_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[2];
    }

    /**
     * @brief Variable for the pion softmax score of the particle.
     * @details The pion softmax score represents the confidence that the
     * network has in the particle being a pion. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a pion.
     * @param p the particle to apply the variable on.
     * @return the pion softmax score of the particle.
     */
    double pion_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[3];
    }

    /**
     * @brief Variable for the proton softmax score of the particle.
     * @details The proton softmax score represents the confidence that the
     * network has in the particle being a proton. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a proton.
     * @param p the particle to apply the variable on.
     * @return the proton softmax score of the particle.
     */
    double proton_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[4];
    }

    /**
     * @brief Variable for the "MIP" softmax score of the particle.
     * @details The "MIP" softmax score is calculated as the sum of the softmax
     * scores for the muon and pion. The score represents the confidence that
     * the network has in the particle being a minimum ionizing particle.
     * @param p the particle to apply the variable on.
     * @return the "MIP" softmax score of the particle.
     */
    double mip_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[2] + p.pid_scores[3];
    }

    /**
     * @brief Variable for the "hadron" softmax score of the particle.
     * @details The "hadron" softmax score is calculated as the sum of the softmax
     * scores for the pion and proton. The score represents the confidence that
     * the network has in the particle being a hadron.
     * @param p the particle to apply the variable on.
     * @return the "hadron" softmax score of the particle.
     */
    double hadron_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[3] + p.pid_scores[4];
    }

    /**
     * @brief Variable for the primary softmax score of the particle.
     * @details The primary softmax score represents the confidence that the
     * network has in the particle being a primary particle. The score is between
     * 0 and 1, with 1 being the most confident that the particle is a primary
     * particle.
     * @param p the particle to apply the variable on.
     * @return the primary softmax score of the particle.
     */
    double primary_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.primary_scores[1];
    }

    /**
     * @brief Variable for the secondary softmax score of the particle.
     * @details The secondary softmax score represents the confidence that the
     * network has in the particle being a secondary particle. The score is between
     * 0 and 1, with 1 being the most confident that the particle is a secondary
     * particle.
     * @param p the particle to apply the variable on.
     * @return the secondary softmax score of the particle.
     */
    double secondary_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.primary_scores[0];
    }
}
#endif // PARTICLE_VARIABLES_H
