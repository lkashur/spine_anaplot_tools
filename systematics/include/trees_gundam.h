/**
 * @file trees_gundam.h
 * @brief Header file for the trees namespace.
 * @details This file contains the header for the trees namespace. The trees
 * namespace contains functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework. Different "copying" actions can
 * be performed on the TTrees, such as a simple copy or adding systematics to
 * the output file based on the selected signal candidates and the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#ifndef TREES_GUNDAM_H
#define TREES_GUNDAM_H
#include <iostream>
#include <sstream>

#include "utilities.h"
#include "configuration.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TClonesArray.h"
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
namespace sys::trees_gundam
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
        /**
         * @brief Create the output subdirectory following the nesting outlined
         * in the configuration file.
         */
        TDirectory * directory = (TDirectory *) output;
        directory = create_directory(directory, table.get_string_field("destination").c_str());
        
        /**
         * @brief Create the output TTree with the name specified in the
         * configuration file.
         */
        TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());

        /**
         * @brief Connect to the input TTree and associated branches.
         * @details Three are N+3 branches in the input TTree, where N is the
         * number of branches in the input TTree of type double. The three
         * other branches (Run, Subrun, Evt) are of type int. This block uses
         * a single array to store the values of the double branches and three
         * separate variables to store the values of the int branches.
         */
        TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
        int run, subrun, event;
        double br[input_tree->GetNbranches()-3];
        for (int i = 0; i < input_tree->GetNbranches()-3; i++)
            input_tree->SetBranchAddress(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        input_tree->SetBranchAddress("Run", &run);
        input_tree->SetBranchAddress("Subrun", &subrun);
        input_tree->SetBranchAddress("Evt", &event);

        /**
         * @brief Create the branches in the output TTree following the same
         * structure as the input TTree.
         * @details The same array (type double) and three variables (type int)
         * as used above with the input TTree are used to create the branches
         * in the output TTree. The branches are created in the same order as
         * the input TTree. This process streamlines the copying of the values
         * from the input TTree to the output TTree.
         */
        for (int i = 0; i < input_tree->GetNbranches()-3; i++)
            output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        output_tree->Branch("Run", &run);
        output_tree->Branch("Subrun", &subrun);
        output_tree->Branch("Evt", &event);

        /**
         * @brief Loop over the input TTree and copy the values to the output
         * TTree.
         */
        for(int i(0); i < input_tree->GetEntries(); ++i)
        {
            input_tree->GetEntry(i);
            output_tree->Fill();
        }

        /**
         * @brief Write the output TTree to the output ROOT file.
         */
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
        /**
         * @brief Create the output subdirectory following the nesting outlined
         * in the configuration file.
         */
        TDirectory * directory = (TDirectory *) output;
        directory = create_directory(directory, table.get_string_field("destination").c_str());
        
        /**
         * @brief Connect to the input TTree and associated branches.
         * @details Three are N+3 branches in the input TTree, where N is the
         * number of branches in the input TTree of type double. The three
         * other branches (Run, Subrun, Evt) are of type int. This block uses
         * a single array to store the values of the double branches and three
         * separate variables to store the values of the int branches. There is
         * one quirk, however, as we would also like to have access to the
         * "nu_id" branch in the input TTree directly.
         */
        TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
	TH1D * _POT = (TH1D *) input->Get("events/mc/POT");
	TH1D *_Livetime = (TH1D*) input->Get("events/mc/Livetime");
        double br[input_tree->GetNbranches()-3];
        double nu_id;
        Int_t run, subrun, event;
        for(int i(0); i < input_tree->GetNbranches()-3; ++i)
            input_tree->SetBranchAddress(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
        input_tree->SetBranchAddress("nu_id", &nu_id);
        input_tree->SetBranchAddress("Run", &run);
        input_tree->SetBranchAddress("Subrun", &subrun);
        input_tree->SetBranchAddress("Evt", &event);

	// Recast variables for GUNDAM input
        Double_t cuttype;
	Double_t issignal;
	Double_t isdata;
        input_tree->SetBranchAddress("CutType", &cuttype);
	input_tree->SetBranchAddress("IsSignal", &issignal);
	input_tree->SetBranchAddress("IsData", &isdata);
	
        /**
         * @brief Create the output TTree with the name specified in the
         * configuration file.
         * @details The output TTree is created with the same branches as the
         * input TTree, plus the Run, Subrun, and Evt branches. The same array
         * (type double) and three variables (type int) as used above with the
         * input TTree are used to create the branches in the output TTree. The
         * branches are created in the same order as the input TTree. This
         * process streamlines the copying of the values from the input TTree to
         * the output TTree.
         */
        //TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());
	TTree * output_tree = new TTree("SelectedEvents_NSigmas", "SelectedEvents_NSigmas");
	TH1D *POT = (TH1D*) _POT->Clone("POT");
	TH1D *Livetime = (TH1D*) _Livetime->Clone("Livetime");
	directory->WriteObject(POT, "POT");
	directory->WriteObject(Livetime, "Livetime");
        for (int i(0); i < input_tree->GetNbranches()-3; ++i)
	{
	  //output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
	  std::string branch = input_tree->GetListOfBranches()->At(i)->GetName();
	    if(strcmp(branch.c_str(), "CutType") && strcmp(branch.c_str(), "IsSignal") && strcmp(branch.c_str(), "IsData"))
	    {
	        output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
	    }
	}
        output_tree->Branch("Run", &run);
        output_tree->Branch("Subrun", &subrun);
        output_tree->Branch("Evt", &event);

	// Recast variables for GUNDAM input
	Int_t _cuttype;
	Int_t _issignal;
	Int_t _isdata;
	output_tree->Branch("CutType", &_cuttype);
	output_tree->Branch("IsSignal", &_issignal);
	output_tree->Branch("IsData", &_isdata);
	
        /**
         * @brief Create the map of selected signal candidates.
         * @details This block creates a map of selected signal candidates. The
         * map is built by looping over the input TTree and storing the index
         * (run, subrun, event, nu_id) of the selected signal candidates in the
         * map. The index is used to match the selected signal candidates with
         * the universe weights for the parent neutrino.
         */
        map_t candidates;
        for(int i(0); i < input_tree->GetEntries(); ++i)
        {
            input_tree->GetEntry(i);
            candidates.insert(std::make_pair<index_t, size_t>(std::make_tuple(run, subrun, event, nu_id), i));
        }

        /**
         * @brief Configure the weight-based systematics.
         * @details This block configures the weight-based systematics. The
         * systematics are split (by type) into separate TTrees, which is
         * enforced by the configuration file. Because we do not wish to loop
         * over the selected signal candidates multiple times, we must store
         * the systematic information in such a way that we can easily
         * accomodate this scheme. The variable "systables" is a vector of
         * vectors of @ref sys::cfg::ConfigurationTables, with an inner vector
         * for each systematic parameter and an outer vector for grouping the
         * systematics by type. The variable "systs" is a map that maps the
         * name of the systematic parameter to its index within the input
         * weights. The variable "systrees" is a map that serves as a container
         * for the output TTrees of each systematic type keyed by the name of
         * the type. The variable "weights" is a map with a key of the index of
         * the systematic parameter and a value of std::vector<double>* that is
         * used to connect the universe weights to the output TTree.
         */
        std::vector<std::vector<sys::cfg::ConfigurationTable>> systables;
	std::map<std::string, TTree *> systrees;
	// Multisim
	std::map<std::string, int64_t> systsMultisim;
	std::map<int64_t, std::vector<double>*> weightsMultisim;
	// Multisigma
	int _sigmas[] = {-1, 1, -2, 2, -3, 3, 0};
	int _sigmas_size = sizeof(_sigmas) / sizeof(_sigmas[0]);
        std::map<std::string, int64_t> systsMultisigma;
	std::map<int64_t, std::vector<double>*> sigmasMultisigma;
	std::map<int64_t, std::vector<double>*> weightsMultisigma;
	TClonesArray *arrMultisigma[config.get_subtables("multisigma").size()];
  
        /**
         * @brief Loop over the systematic types in the configuration file.
         * @details This block loops over the systematic types in the
         * configuration file. The "table" field in the configuration file
         * specifies the name of the table lists in the configuration file that
         * contain the exact definition of the systematics. Principally, this
         * loop is used to load and configure the systematics of each type in
         * sequential order.
         */
        for(std::string & s : table.get_string_vector("table"))
        {
            systables.push_back(config.get_subtables(s));
            systrees[s] = new TTree((s+"Tree").c_str(), (s+"Tree").c_str());
	    //systrees["multisim"] = new TTree((s+"Tree").c_str(), (s+"Tree").c_str());

            /**
             * @brief Loop over the systematics of this type.
             * @details This block loops over systematics belonging to the same
             * type as defined by the configuration file. The loop "flattens"
             * the systematics into a single set of maps for use below in the
             * loop over the input CAF files. The TTree of each systematic type
             * is extended with a vector of doubles for each systematic
             * parameter belonging to the type.
             */
            for(sys::cfg::ConfigurationTable & t : systables.back())
            {
	      
	        // Multisim
	        if(!strcmp(s.c_str(), "multisim"))
		{
		    systsMultisim.insert(std::make_pair<std::string, int64_t>(t.get_string_field("name"), t.get_int_field("index")));
		    weightsMultisim.insert(std::make_pair<int64_t, std::vector<double>*>(t.get_int_field("index"), new std::vector<double>));
		    systrees[s]->Branch(t.get_string_field("name").c_str(), &weightsMultisim[t.get_int_field("index")]);
		}
		// Multisigma
		else if(!strcmp(s.c_str(), "multisigma"))
		{
		    systsMultisigma.insert(std::make_pair<std::string, int64_t>(t.get_string_field("name_short"), t.get_int_field("index")));
		    sigmasMultisigma.insert(std::make_pair<int64_t, std::vector<double>*>(t.get_int_field("index"), new std::vector<double>));
		    weightsMultisigma.insert(std::make_pair<int64_t, std::vector<double>*>(t.get_int_field("index"), new std::vector<double>));
		    //systrees[s]->Branch(t.get_string_field("name").c_str(), &weightsMultisigma[t.get_int_field("index")]);
		}  
            }
        }

	int sysIdx1(0);
        for(auto & [key, value] : systsMultisigma)
	{  
	    arrMultisigma[sysIdx1] = new TClonesArray("TGraph", 1);
	    //systrees["multisigma"]->Branch((keyShort).c_str(), &arrMultisigma[sysIdx1], 32000, -1);
	    output_tree->Branch(key.c_str(), &arrMultisigma[sysIdx1], 32000, -1);
	    sysIdx1++;
	}

        /**
         * @brief Load the input CAF files.
         * @details This block loads the input CAF files. The input CAF files
         * are specified in a file list that is read from the configuration
         * file. The file list is read line-by-line and the input CAF files are
         * stored in a vector.
         */
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
			_cuttype = (Int_t) cuttype;
			_issignal = (Int_t) issignal;
			_isdata = (Int_t) isdata;
                        //output_tree->Fill();
                        
                        /**
                         * @brief Store the universe weights in the output TTree.
                         * @details This block stores the universe weights in the
                         * output TTree for each of the configured systematics.  
                         */
                        for(auto & [key, value] : systsMultisim)
                        {
                            weightsMultisim[value]->clear();
                            for(size_t u(0); u < nu.wgt[value].univ.size(); ++u)
                                weightsMultisim[value]->push_back(nu.wgt[value].univ[u]);
                        } // End of loop over the configured systematics.
			
			int sysIdx2(0);
			for(auto & [key, value] : systsMultisigma)
			{
			    sigmasMultisigma[value]->clear();
			    weightsMultisigma[value]->clear();
			    
			    // Sigmas
			    for(int i(0); i < _sigmas_size; i++)
			    {
			        sigmasMultisigma[value]->push_back(_sigmas[i]);
			    }
			    
			    // Weights
			    int numWgts(0);
			    for(size_t u(0); u < nu.wgt[value].univ.size(); ++u)
			    {
			        weightsMultisigma[value]->push_back(nu.wgt[value].univ[u]);
				numWgts++;
			    }
			    weightsMultisigma[value]->push_back(1);

			    // Fill TGraph
			    double* sigmasArr = sigmasMultisigma[value]->data();
                            double* weightsArr = weightsMultisigma[value]->data();
                            new( (*arrMultisigma[sysIdx2])[0]) TGraph(numWgts+1,sigmasArr,weightsArr);
                            sysIdx2++;

			} // End of loop over the configured systematics.
			output_tree->Fill();
                        for(auto & [key, value] : systrees)
                            value->Fill();
                    }
                } // End of loop over the neutrino interactions in the input CAF file.
            } // End of loop over the events in the input CAF file.
            caf->Close();
            nprocessed++;
        } // End of loop over the input CAF files.
        input->Close();
        //directory->WriteObject(output_tree, table.get_string_field("name").c_str());
	directory->WriteObject(output_tree, "SelectedEvents_NSigmas");
        for(auto & [key, value] : systrees)
            directory->WriteObject(value, (key+"Tree").c_str());
    }   
}
#endif
