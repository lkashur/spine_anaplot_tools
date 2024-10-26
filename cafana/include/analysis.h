/**
 * @file analysis.h
 * @brief Header file for the Analysis class, a class designed to streamline
 * the running of multiple samples through CAFAna.
 * @details The Analysis class operates under the principle that the full
 * dataset in an analysis is comprised of multiple samples, each of which
 * represents a different component of the signal or background. Under this
 * paradigm, an analysis consists of running the same set of variables/cuts
 * on each sample, and then combining the results at the end. The Analysis
 * class is designed to facilitate this process by allowing the user to
 * specify the variables and cuts to be applied to each sample, along with the
 * full list of samples to be run. The class then handles the running of the
 * samples and the storing of the results.
 * @author mueller@fnal.gov
 */
#ifndef ANALYSIS_H
#define ANALYSIS_H
#include <vector>
#include <string>

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "TDirectory.h"
#include "TFile.h"

/**
 * @namespace ana
 * @brief Namespace for the Analysis class and related functions.
 * @details The ana namespace contains the Analysis class, which is designed
 * to streamline the running of multiple samples through CAFAna. The namespace
 * is also used to organize analysis-related functions and variables within
 * CAFAna.
 */
namespace ana
{
    /**
     * @class Analysis
     * @brief Class designed to streamline the running of multiple samples
     * through CAFAna.
     * @details The Analysis class operates under the principle that the full
     * dataset in an analysis is comprised of multiple samples, each of which
     * represents a different component of the signal or background. Under this
     * paradigm, an analysis consists of running the same set of variables/cuts
     * on each sample, and then combining the results at the end. The Analysis
     * class is designed to facilitate this process by allowing the user to
     * specify the variables and cuts to be applied to each sample, along with
     * the full list of samples to be run. The class then handles the running of
     * the samples and the storing of the results.
     */
    class Analysis
    {
        public:
            Analysis(std::string name);
            void AddLoader(std::string name, ana::SpectrumLoader * loader);
            void AddVars(std::vector<std::string> names, std::vector<ana::SpillMultiVar> & vars);
            void Go();
        private:
            std::string name;
            std::vector<std::string> lnames;
            std::vector<ana::SpectrumLoader*> loaders;
            std::vector<std::string> names;
            std::vector<ana::SpillMultiVar> vars;
    };

    /**
     * @brief Default constructor for the Analysis class.
     * @details This constructor initializes an instance of the Analysis class
     * with a default name. 
     * @param name The name of the Analysis class. This name is used in the
     * creation of the output ROOT file.
     * @return A new instance of the Analysis class.
     */
    Analysis::Analysis(std::string name)
    {
        this->name = name;
    }

    /**
     * @brief Add a SpectrumLoader to the Analysis class.
     * @details This function allows the user to add a SpectrumLoader to the
     * Analysis class. The SpectrumLoader represents a sample in the analysis,
     * and is used to load the data from the ROOT file and apply the cuts and
     * variables to the data.
     * @param name The name of the SpectrumLoader to be added to the Analysis
     * class.
     * @param loader The SpectrumLoader to be added to the Analysis class.
     * @return void
     */
    void Analysis::AddLoader(std::string name, ana::SpectrumLoader * loader)
    {
        lnames.push_back(name);
        loaders.push_back(loader);
    }

    /**
     * @brief Add a set of variables to the Analysis class.
     * @details This function allows the user to add a set of variables to the
     * Analysis class. The variables are represented by a vector of SpillMultiVar
     * objects, which contain the names of the variables and the cuts to be
     * applied to the data.
     * @param names A vector of strings containing the names of the variables to
     * be added to the Analysis class.
     * @param vars A vector of SpillMultiVar objects containing the variables and
     * cuts to be added to the Analysis class.
     * @return void
     */
    void Analysis::AddVars(std::vector<std::string> names, std::vector<ana::SpillMultiVar> & vars)
    {
        this->names.insert(this->names.end(), names.begin(), names.end());
        this->vars.insert(this->vars.end(), vars.begin(), vars.end());
    }

    /**
     * @brief Run the analysis on the specified samples.
     * @details This function runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample, applies the cuts and variables to the data, and stores the
     * results in a TFile.
     * @return void
     */
    void Analysis::Go()
    {
        TFile * f = new TFile(std::string(name + ".root").c_str(), "RECREATE");
        TDirectory * dir = f->mkdir("events");
        dir->cd();

        for(size_t i(0); i < loaders.size(); ++i)
        {
            TDirectory * subdir = dir->mkdir(lnames[i].c_str());
            subdir->cd();
            ana::Tree tree = ana::Tree("selectedNu", names, *loaders[i], vars, ana::kNoSpillCut, true);
            loaders[i]->Go();
            tree.SaveTo(subdir);
            dir->cd();
        }
        f->Close();
    }
}

#endif // ANALYSIS_H