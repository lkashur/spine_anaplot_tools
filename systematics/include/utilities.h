/**
 * @file utilities.h
 * @brief Header and implementation of functions generally useful for the
 * systematics framework within SPINE.
 * @details This file contains the implementation of functions that are
 * generally useful for the systematics framework within SPINE. These functions
 * are used to read in the selected signal candidates.
 * @author mueller@fnal.gov
 */
#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include <map>
#include <string>
#include <tuple>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

typedef std::tuple<Double_t, Double_t, Double_t, Double_t> index_t;
typedef std::map<index_t, size_t> map_t;

/**
 * @brief Read the selected signal candidates from the TTree.
 * @details This function reads the selected signal candidates from the TTree
 * and stores them in a map. The map is indexed by the run, subrun, event, and
 * nu_id of the candidate. The value of the map is the entry number of the
 * candidate in the TTree. This map can be used later to match back to the
 * orignal candidates when writing the output TTree.
 * @param candidates The map that will store the selected signal candidates.
 * @param reader The TTreeReader that is used to read the TTree.
 * @return void
 */
void read_candidates(map_t & candidates, TTreeReader & reader)
{
    TTreeReaderValue<Int_t> run(reader, "Run");
    TTreeReaderValue<Int_t> subrun(reader, "Subrun");
    TTreeReaderValue<Int_t> event(reader, "Evt");
    TTreeReaderValue<Double_t> nu_id(reader, "nu_id");

    while(reader.Next())
        candidates.insert(std::make_pair<index_t, size_t>(std::make_tuple(*run, *subrun, *event, *nu_id), reader.GetCurrentEntry()));
}

/**
 * @brief Create a nested directory structure in the output ROOT file.
 * @details This function creates a nested directory structure in the output
 * ROOT file. The directory structure is parsed from the desired output
 * destination in the configuration file.
 * @param parent The parent directory in the output ROOT file.
 * @param path The desired output destination.
 * @return directory The final directory in the output ROOT file.
 */
TDirectory * create_directory(TDirectory * parent, std::string path)
{
    TDirectory * directory = parent;
    size_t pos(path.find_last_of("/"));
    if(pos != std::string::npos)
    {
        std::string subdirs(path.substr(0, pos));
        std::string subdir;
        while((pos = subdirs.find_first_of("/")) != std::string::npos)
        {
            subdir = subdirs.substr(0, pos);
            if(directory->Get(subdir.c_str()) == 0)
                directory = directory->mkdir(subdir.c_str());
            else
                directory = directory->GetDirectory(subdir.c_str());
            subdirs = subdirs.substr(pos+1);
        }
        if(directory->Get(subdirs.c_str()) == 0)
            directory = directory->mkdir(subdirs.c_str());
        else
            directory = directory->GetDirectory(subdirs.c_str());
    }
    return directory;
}
#endif // UTILITIES_H