/**
 * @file configuration.cc
 * @brief Implementation of the ConfigurationTable class and related functions.
 * @details This file contains the implementation of the ConfigurationTable
 * class and related functions. The ConfigurationTable class is used to read
 * and interface with the TOML-based configuration file. The configuration file
 * consolidates all the parameters that may be analysis-specific or that may
 * change frequently. The main backend employed in this project is the toml++
 * library.
 * @author mueller@fnal.gov
 */
#include <string>

#include "configuration.h"
#include "toml++/toml.h"

namespace sys::cfg
{
    // Configuration class constructor that takes a TOML configuration table
    // and stores it as a class member (for use with sub-tables of the
    // configuration file).
    ConfigurationTable::ConfigurationTable(const toml::table & config) : config(config) {}

    // Set the configuration file by reading the TOML configuration table and
    // validating the configuration file.
    void ConfigurationTable::set_config(const std::string & path)
    {
        try
        {
            config = toml::parse_file(path);
        }
        catch(const std::exception & e)
        {
            throw ConfigurationError(e.what());
        }
        validate();
    }

    // Check that the requested field is present in the configuration file.
    void ConfigurationTable::check_field(const std::string & field)
    {
        std::optional<std::string> value(config.at_path(field).value<std::string>());
        if(!value)
            throw ConfigurationError("Field " + field + " not found in the configuration file.");
    }

    // Check that the requested field is present in the configuration file.
    bool ConfigurationTable::has_field(const std::string & field)
    {
        return config.contains(field);
    }

    // Retrieve the requested string field from the configuration table.
    std::string ConfigurationTable::get_string_field(const std::string & field)
    {
        std::optional<std::string> value(config.at_path(field).value<std::string>());
        if(!value)
            throw ConfigurationError("Field " + field + " (string) not found in the configuration file.");
        return *value;
    }

    // Retrieve the requested vector of strings from the configuration table.
    std::vector<std::string> ConfigurationTable::get_string_vector(const std::string & field)
    {
        std::vector<std::string> values;
        toml::array * elements = config.at_path(field).as_array();
        for(auto & e : *elements)
            values.push_back(*e.value<std::string>());
        return values;
    }

    // Retrieve the requested integer field from the configuration table.
    int64_t ConfigurationTable::get_int_field(const std::string & field)
    {
        std::optional<int64_t> value(config.at_path(field).value<int64_t>());
        if(!value)
            throw ConfigurationError("Field " + field + " (int) not found in the configuration file.");
        return *value;
    }

    // Retrieve the requested double field from the configuration table.
    double ConfigurationTable::get_double_field(const std::string & field)
    {
        std::optional<double> value(config.at_path(field).value<double>());
        if(!value)
            throw ConfigurationError("Field " + field + " (double) not found in the configuration file.");
        return *value;
    }

    // Retrieve the requested vector of doubles from the configuration table.
    std::vector<double> ConfigurationTable::get_double_vector(const std::string & field)
    {
        std::vector<double> values;
        toml::array * elements = config.at_path(field).as_array();
        for(auto & e : *elements)
            values.push_back(*e.value<double>());
        return values;
    }

    // Validate the configuration file by checking that all the requisite
    // fields are present.
    void ConfigurationTable::validate()
    {
        check_field("input.path");
        check_field("output.path");
    }

    // Get a list of all subtables matching the requested table name.
    std::vector<ConfigurationTable> ConfigurationTable::get_subtables(const std::string & table)
    {
        std::vector<ConfigurationTable> tables;
        if(config.find(table) == config.end())
            throw ConfigurationError("Table " + table + " not found in the configuration file.");
        toml::array * elements = config[table].as_array();
        for(auto & e : *elements)
            tables.push_back(ConfigurationTable(*e.as_table()));
        return tables;
    }
} // namespace systematics::configuration