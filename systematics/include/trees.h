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

#include "utilities.h"
#include "configuration.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

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
     * @brief Type definitions for the selected signal candidates and the
     * universe weights.
     */
    typedef std::tuple<Double_t, Double_t, Double_t, Double_t> index_t;
    typedef std::map<index_t, size_t> map_t;

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
    void copy_tree(sys::cfg::ConfigurationTable & table, TFile * output, TFile * input)
    {
        // Create the output subdirectory.
        TDirectory * directory = (TDirectory *) output;
        directory = create_directory(directory, table.get_string_field("destination").c_str());
        
        // Create the output TTree.
        TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());

        // Connect to the input TTree.
        TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
        int run, subrun, event;
        double br[input_tree->GetNbranches()-3];
        for (int i = 0; i < input_tree->GetNbranches()-3; i++)
            input_tree->SetBranchAddress(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        input_tree->SetBranchAddress("Run", &run);
        input_tree->SetBranchAddress("Subrun", &subrun);
        input_tree->SetBranchAddress("Evt", &event);

        // Create the output TTree branches.
        for (int i = 0; i < input_tree->GetNbranches()-3; i++)
            output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        output_tree->Branch("Run", &run);
        output_tree->Branch("Subrun", &subrun);
        output_tree->Branch("Evt", &event);

        // Loop over the input TTree and copy the values to the output TTree.
        for(int i(0); i < input_tree->GetEntries(); ++i)
        {
            input_tree->GetEntry(i);
            output_tree->Fill();
        }

        // Write the output TTree to the output file.
        directory->WriteObject(output_tree, table.get_string_field("name").c_str());
    }

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
    void copy_with_weight_systematics(sys::cfg::ConfigurationTable & config, sys::cfg::ConfigurationTable & table, TFile * output, TFile * input)
    {
        // Create the output subdirectory.
        TDirectory * directory = (TDirectory *) output;
        directory = create_directory(directory, table.get_string_field("destination").c_str());
        
        // Connect to the input TTree.
        TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
        double br[input_tree->GetNbranches()-3];
        double nu_id;
        Int_t run, subrun, event;
        for(int i(0); i < input_tree->GetNbranches()-3; ++i)
            input_tree->SetBranchAddress(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        input_tree->SetBranchAddress("nu_id", &nu_id);
        input_tree->SetBranchAddress("Run", &run);
        input_tree->SetBranchAddress("Subrun", &subrun);
        input_tree->SetBranchAddress("Evt", &event);

        // Create and populate a map of the selected signal candidates.
        map_t candidates;
        for(int i(0); i < input_tree->GetEntries(); ++i)
        {
            input_tree->GetEntry(i);
            candidates.insert(std::make_pair<index_t, size_t>(std::make_tuple(run, subrun, event, nu_id), i));
        }

        // Create the output TTree for the selected signal candidates.
        TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());
        for (int i(0); i < input_tree->GetNbranches()-3; ++i)
            output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        
        output_tree->Branch("Run", &run);
        output_tree->Branch("Subrun", &subrun);
        output_tree->Branch("Evt", &event);

        // Configure the weight-based systematics.
        std::vector<std::vector<sys::cfg::ConfigurationTable>> systables;
        std::map<std::string, int64_t> systs;
        std::map<std::string, TTree *> systrees;
        std::map<int64_t, std::vector<double>*> weights;
        for(std::string & s : table.get_string_vector("table"))
        {
            systables.push_back(config.get_subtables(s));
            systrees[s] = new TTree((s+"Tree").c_str(), (s+"Tree").c_str());
            for(sys::cfg::ConfigurationTable & t : systables.back())
            {
                systs.insert(std::make_pair<std::string, int64_t>(t.get_string_field("name"), t.get_int_field("index")));
                weights.insert(std::make_pair<int64_t, std::vector<double>*>(t.get_int_field("index"), new std::vector<double>));
                systrees[s]->Branch(t.get_string_field("name").c_str(), &weights[t.get_int_field("index")]);
            }
        }

        // Read in the input list of CAF files.
        std::vector<std::string> input_files;
        std::ifstream file_list(config.get_string_field("input.caflist"));
        std::string line;
        while(std::getline(file_list, line))
            input_files.push_back(line);
        file_list.close();

        /**
         * @brief Loop over the input CAF files.
         * @details This block loops over the input CAF files. This loop begins
         * the process of matching the selected signal candidates with the universe
         * weights for parent neutrino. 
         */
        size_t nprocessed(0);
        for(std::string input_file : input_files)
        {
            if(nprocessed % 100 == 0)
                std::cout << "Processed " << nprocessed << " files." << std::endl;
            /**
             * @brief Open and validate the input CAF file.
             * @details This block opens the input CAF file and checks that the
             * files is not a zombie (corrupted) and that it contains the TTree
             * "recTree". If the file is a zombie or does not contain the TTree,
             * the code prints an error message and continues to the next file.
             * Practically, skipping a file means that signal events from that file
             * will not be included in the output ROOT file, though a situation
             * where this occurs is not expected.
             */
            TFile * caf = TFile::Open(input_file.c_str(), "READ");
            if(caf->IsZombie() || !caf->GetListOfKeys()->Contains("recTree"))
            {
                std::cerr << "Error: File" << input_file << " does not exist." << std::endl;
                continue;
            }

            /**
             * @brief Connect to the input CAF file fields.
             * @details This block connects to a minimal set of the input CAF file
             * fields. We need the run, subrun, event, and the true interaction
             * (parent neutrino) information to match the selected signal
             * candidates and retrieve the universe weights.
             */
            TTreeReader reader("recTree", caf);
            TTreeReaderValue<uint32_t> rrun(reader, "rec.hdr.run");
            TTreeReaderValue<uint32_t> rsubrun(reader, "rec.hdr.subrun");
            TTreeReaderValue<uint32_t> revt(reader, "rec.hdr.evt");
            TTreeReaderArray<caf::SRTrueInteraction> mc(reader, "rec.mc.nu");

            /**
             * @brief Loop over the events in the input CAF file.
             * @details This block loops over the events in the input CAF file. At
             * each event, the code checks if a neutrino interaction from the input
             * CAF file matches a selected signal candidate. If a match is found,
             * the code retrieves the universe weights for the parent neutrino and
             * stores them in the output TTree.
             */
            while(reader.Next())
            {
                for(const caf::SRTrueInteraction & nu : mc)
                {
                    index_t index(*rrun, *rsubrun, *revt, nu.index);
                    if(candidates.find(index) != candidates.end())
                    {
                        /**
                         * @brief Retrieve the selected signal candidate and copy
                         * the values to the output TTree.
                         * @details This block retrieves the selected signal
                         * candidate that has been matched with the parent neutrino
                         * and copies the values to the output TTree.
                         */
                        input_tree->GetEntry(candidates[index]);
                        run = *rrun;
                        subrun = *rsubrun;
                        event = *revt;
                        output_tree->Fill();
                        
                        /**
                         * @brief Store the universe weights in the output TTree.
                         * @details This block stores the universe weights in the
                         * output TTree for each of the configured systematics.  
                         */
                        for(auto & [key, value] : systs)
                        {
                            weights[value]->clear();
                            for(size_t u(0); u < nu.wgt[value].univ.size(); ++u)
                                weights[value]->push_back(nu.wgt[value].univ[u]);
                        } // End of loop over the configured systematics.
                        for(auto & [key, value] : systrees)
                            value->Fill();
                    }
                } // End of loop over the neutrino interactions in the input CAF file.
            } // End of loop over the events in the input CAF file.
            caf->Close();
            nprocessed++;
        } // End of loop over the input CAF files.
        input->Close();
        directory->WriteObject(output_tree, table.get_string_field("name").c_str());
        for(auto & [key, value] : systrees)
            directory->WriteObject(value, (key+"Tree").c_str());
    }   
}

#endif