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