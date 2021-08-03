#include "file.h"
#include "helpers.h"
#include "interpro.h"
#include "protein.h"

namespace mipfinder::interpro
{

    //Parse the Interpro entry list (available from https://www.ebi.ac.uk/interpro/download/.
    //The format of the data should be a tab-separated file consisting of interpro accession number,
    //entry type and entry name. E.g. "IPR000126	Active_site	Serine proteases, V8 family, serine active site"
    //The first line of the file should be the three headers for the columns. If the file
    //does not correspond to the given, the behaviour is undefined.
    Entries parseEntryList(const std::filesystem::path& interpro_entry_list)
    {
        mipfinder::interpro::Entries results;
        auto stream = mipfinder::file::open(interpro_entry_list);
        std::string line;
        std::getline(stream, line); /* Skips the header in database file */
        while (std::getline(stream, line)) {
            const auto tokens = mipfinder::tokenise(line, '\t');

            /* Correctly formatted file has three columns of data */
            if (tokens.size() != 3) {
                continue;
            }

            mipfinder::interpro::Type type;
            const auto& domain_type = tokens[1];
            if (domain_type == "Active_site") {
                type = mipfinder::interpro::Type::active_site;
            }
            else if (domain_type == "Binding_site") {
                type = mipfinder::interpro::Type::binding_site;
            }
            else if (domain_type == "Conserved_site") {
                type = mipfinder::interpro::Type::conserved_site;
            }
            else if (domain_type == "Domain") {
                type = mipfinder::interpro::Type::domain_type;
            }
            else if (domain_type == "Family") {
                type = mipfinder::interpro::Type::family;
            }
            else if (domain_type == "Homologous_superfamily") {
                type = mipfinder::interpro::Type::homologous_superfamily;
            }
            else if (domain_type == "PTM") {
                type = mipfinder::interpro::Type::ptm;
            }
            else if (domain_type == "Repeat") {
                type = mipfinder::interpro::Type::repeat;
            }
            else {
                type = mipfinder::interpro::Type::unknown;
            }

            const auto& interpro_accession = tokens[0];
            const auto& entry_description = tokens[2];
            results[interpro_accession] = mipfinder::interpro::Entry{.description = entry_description, .type = type};
        }
        return results;
    }

    ProteinDomains parseProteinDomainList(const std::filesystem::path& uniprot_to_interpro_table)
    {
        auto stream = mipfinder::file::open(uniprot_to_interpro_table);
        std::string line;
        mipfinder::interpro::ProteinDomains results;
        while (std::getline(stream, line)) {
            const auto tokens = mipfinder::tokenise(line, '\t');

            /* Correctly formatted file has three columns of data */
            if (tokens.size() != 3) {
                continue;
            }

            std::string uniprot_id = tokens[0];
            std::string sequence_version = tokens[1];
            std::string full_protein_id = uniprot_id + sequence_version;
            constexpr auto domain_entry_delimiter = mipfinder::Protein::id_delimiter;
            std::vector<std::string> interpro_identifiers = mipfinder::tokenise(tokens[2], domain_entry_delimiter);
            for (const auto& identifier : interpro_identifiers) {
                results[full_protein_id].insert(identifier);
            }
        }
        return results;
    }
}