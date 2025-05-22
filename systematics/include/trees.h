/**
 * @file trees.h
 * @brief Header file for the trees namespace.
 * @details This file contains the header for the trees namespace. The trees
 * namespace contains functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework. Different "copying" actions can
 * be performed on the TTrees, such as a simple copy or adding systematics to
 * the output file based on the selected signal candidates and the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#ifndef TREES_H
#define TREES_H
#include <iostream>

#include "detsys.h"
#include "configuration.h"

#include "TFile.h"

/**
 * @namespace sys::trees
 * @brief Namespace for functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework.
 * @details This namespace contains functions that read and interface with the
 * TTrees produced by the CAFAna analysis framework. Different "copying" 
 * actions can be performed on the TTrees, such as a simple copy or the
 * addition of systematics to the output file based on the selected signal
 * candidates and the configured systematics.
 */
namespace sys::trees
{
    /**
     * @brief Type definitions for systematic indexing (variable name and
     * index).
     */
    typedef std::pair<std::string, int64_t> syst_t;

    /**
     * @brief Simple hash function for a set of five variables.
     * @details This function takes five variables (run, subrun, event,
     * nu_id, and nu_energy) and uses them to create a unique index value. This
     * is done by bit-packing the variables into a single 64-bit unsigned
     * integer.
     * @param run The run number.
     * @param subrun The subrun number.
     * @param event The event number.
     * @param nu_id The neutrino ID.
     * @param nu_energy The neutrino energy (default is 0).
     * @return The hash value.
     */
    size_t hash(uint64_t run, uint64_t subrun, uint64_t event, uint64_t nu_id, float nu_energy=0);

    /**
     * @brief Copy the input TTree to the output TTree.
     * @details This function copies the input TTree to the output TTree. The
     * function loops over the input TTree and copies the values of the branches
     * to the output TTree. The output TTree is created with the same branches
     * as the input TTree.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_tree(sys::cfg::ConfigurationTable & table, TFile * output, TFile * input);

    /**
     * @brief Add reweightable systematics to the output TTree.
     * @details This function adds reweightable systematics to the output
     * TTree. The function loops over the input TTree to build a map for the
     * selected signal candidates to their index in the input TTree. The
     * function then loops over the neutrinos in the CAF input files and
     * populates the output TTree with the selected signal candidates and the
     * universe weights for matched neutrinos.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_with_weight_systematics(sys::cfg::ConfigurationTable & config, sys::cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc);
}
#endif