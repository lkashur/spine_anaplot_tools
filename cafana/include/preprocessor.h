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
 * Preprocessor wrapper for looping over reco interactions and broadcasting a
 * function over the matching true interactions. The SpillMultiVar accepts a
 * vector as a result of some function running over the top-level
 * StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each true interaction that
 * is matched to by a reco interaction passing the cut SEL.
*/
#define VARDLP_RECO_TRUE(NAME,VAR,SEL)                       \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr) \
    {							                             \
        std::vector<double> var;				             \
        for(auto const& i : sr->dlp)			             \
        {							                         \
            if(SEL(i) && i.match.size() > 0)		         \
	        var.push_back(VAR(sr->dlp_true[i.match[0]]));    \
        }							                         \
        return var;						                     \
  })

/**
 * Preprocessor wrapper for looping over reco interactions and broadcasting a
 * function over the reco interactions. The SpillMultiVar accepts a vector as
 * a result of some function running over the top-level StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each reco interaction
 * passing the cut SEL.
*/
#define VARDLP_RECO_RECO(NAME,VAR,SEL)                       \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr) \
    {							                             \
        std::vector<double> var;				             \
        for(auto const& i : sr->dlp)			             \
        {							                         \
            if(SEL(i) && i.match.size() > 0)		         \
	        var.push_back(VAR(i));				             \
        }							                         \
        return var;						                     \
    })

/**
 * Preprocessor wrapper for looping over true interactions and broadcasting a
 * function over the true interactions. The SpillMultiVar accepts a vector as
 * a result of some function running over the top-level StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each true interaction
 * passing the cut SEL.
*/
#define VARDLP_TRUE_TRUE(NAME,VAR,SEL)                       \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr) \
    {							                             \
        std::vector<double> var;				             \
        for(auto const& i : sr->dlp_true)			         \
        {							                         \
            if(SEL(i) && i.match.size() > 0)		         \
            var.push_back(VAR(i));                           \
        }							                         \
        return var;						                     \
    })

/**
 * Preprocessor wrapper for looping over true interactions and broadcasting a
 * function over the matching reco interactions. The SpillMultiVar accepts a
 * vector as a result of some function running over the top-level StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each reco interaction that
 * is matched to by a true interaction passing the cut SEL.
*/
#define VARDLP_TRUE_RECO(NAME,VAR,SEL)                       \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr) \
    {							                             \
        std::vector<double> var;				             \
        for(auto const& i : sr->dlp_true)			         \
        {							                         \
            if(SEL(i) && i.match.size() > 0)		         \
            var.push_back(VAR(sr->dlp[i.match[0]]));         \
        }							                         \
        return var;						                     \
    })
#endif // PREPROCESSOR_H