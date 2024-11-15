/**
 * @file preprocessor.h
 * @brief Header file for preprocessor macros that streamline variable
 * declarations.
 * @author mueller@fnal.gov
*/
#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H
#include "sbnana/CAFAna/Core/MultiVar.h"

/**
 * @brief Preprocessor wrapper for looping over reco interactions and
 * broadcasting a function over the reco interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the reco interactions and returns a vector of the results of the
 * function.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @param CAT function broadcast over the true interactions matched to by the
 * reco interactions to select only reco interactions belonging to a certain
 * truth category.
 * @return a vector with the result of VAR called on each reco interaction
 * passing the cut SEL and which belongs to the truth category established by
 * CAT.
 */
#define SPINEVAR_RR(VAR,SEL,CAT)                                                            \
    [](const caf::SRSpillProxy* sr)                                                         \
    {                                                                                       \
        std::vector<double> var;                                                            \
        bool is_mc(sr->ndlp_true != 0);                                                     \
        for(auto const& i : sr->dlp)                                                        \
        {                                                                                   \
            if(SEL(i) && ((i.match_ids.size() > 0 && CAT(sr->dlp_true[i.match_ids[0]])) || !is_mc)) \
                var.push_back(VAR(i));                                                      \
        }                                                                                   \
        return var;                                                                         \
    }

 /**
 * @brief Preprocessor wrapper for looping over reco interactions and
 * broadcasting a function over the matching true interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the reco interactions and returns a vector of the results of the
 * function for each true interaction.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @param CAT function broadcast over the true interactions matched to by the
 * reco interactions to select only reco interactions belonging to a certain
 * truth category.
 * @return a vector with the result of VAR called on each true interaction that
 * is matched to by a reco interaction passing the cut SEL and which belongs to
 * the truth category established by CAT.
 */
#define SPINEVAR_RT(VAR,SEL,CAT)                                                            \
    [](const caf::SRSpillProxy* sr)                                                         \
    {                                                                                       \
        std::vector<double> var;                                                            \
        bool is_mc(sr->ndlp_true != 0);                                                     \
        for(auto const& i : sr->dlp)                                                        \
        {                                                                                   \
            if(SEL(i) && ((i.match_ids.size() > 0 && CAT(sr->dlp_true[i.match_ids[0]])) || !is_mc)) \
                var.push_back(i.match_ids.size() > 0 ? VAR(sr->dlp_true[i.match_ids[0]]) : -1.0);   \
        }                                                                                   \
        return var;                                                                         \
 }

/**
 * @brief Preprocessor wrapper for looping over true interactions and
 * broadcasting a function over the matching reco interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the true interactions and returns a vector of the results of the
 * function for each reco interaction.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each reco interaction that
 * is matched to by a true interaction passing the cut SEL.
 */
#define SPINEVAR_TR(VAR,SEL)                         \
    [](const caf::SRSpillProxy* sr)                  \
    {							                     \
        std::vector<double> var;			         \
        for(auto const& i : sr->dlp_true)	         \
        {                                            \
            if(SEL(i) && i.match_ids.size() > 0)         \
            var.push_back(VAR(sr->dlp[i.match_ids[0]])); \
        }                                            \
        return var;                                  \
    }

/**
 * @brief Preprocessor wrapper for looping over true interactions and
 * broadcasting a function over the true interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the true interactions and returns a vector of the results of the
 * function.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each true interaction
 * passing the cut SEL.
 */
#define SPINEVAR_TT(VAR,SEL)                 \
    [](const caf::SRSpillProxy* sr)          \
    {                                        \
        std::vector<double> var;             \
        for(auto const& i : sr->dlp_true)    \
        {                                    \
            if(SEL(i) && i.match_ids.size() > 0) \
                var.push_back(VAR(i));       \
        }                                    \
        return var;                          \
    }
#endif // PREPROCESSOR_H