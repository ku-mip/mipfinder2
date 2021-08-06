#include "helpers.h"
#include "interpro.h"
#include "protein.h"

namespace mipfinder::interpro
{

    //Parse the Interpro entry list (available from https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list [accessed 06/08/21]).
    //The format of the supplied interpro_entry_list should be a tab-separated file consisting of interpro accession number,
    //entry type and entry name. E.g. "IPR000126	Active_site	Serine proteases, V8 family, serine active site"
    //The first line of the file should be the three headers for the columns. If the file
    //does not correspond to the given, the behaviour is undefined.
    //Returns a list of InterPro entries in a non-determinate order.
    Entries parseEntryList(const std::filesystem::path& interpro_entry_list)
    {
        std::ifstream file{ interpro_entry_list };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " + interpro_entry_list.string() + ", aborting...");
        }

        mipfinder::interpro::Entries results;
        std::string line;
        std::getline(file, line); /* Skips the header in database file */
        while (std::getline(file, line)) {
            const auto tokens = mipfinder::tokenise(line, '\t');

            /* Correctly formatted file has three columns of data */
            if (tokens.size() != 3) {
                continue;
            }

            mipfinder::interpro::Entry::Type type;
            const auto& domain_type = tokens[1];
            if (domain_type == "Active_site") {
                type = mipfinder::interpro::Entry::Type::active_site;
            }
            else if (domain_type == "Binding_site") {
                type = mipfinder::interpro::Entry::Type::binding_site;
            }
            else if (domain_type == "Conserved_site") {
                type = mipfinder::interpro::Entry::Type::conserved_site;
            }
            else if (domain_type == "Domain") {
                type = mipfinder::interpro::Entry::Type::domain_type;
            }
            else if (domain_type == "Family") {
                type = mipfinder::interpro::Entry::Type::family;
            }
            else if (domain_type == "Homologous_superfamily") {
                type = mipfinder::interpro::Entry::Type::homologous_superfamily;
            }
            else if (domain_type == "PTM") {
                type = mipfinder::interpro::Entry::Type::ptm;
            }
            else if (domain_type == "Repeat") {
                type = mipfinder::interpro::Entry::Type::repeat;
            }
            else {
                type = mipfinder::interpro::Entry::Type::unknown;
            }

            const auto& interpro_accession = tokens[0];
            const auto& entry_name = tokens[2];
            results.emplace_back(mipfinder::interpro::Entry{.name = entry_name, .accession = interpro_accession, .type = type});
        }
        return results;
    }

    ProteinDomains parseProteinDomainList(const std::filesystem::path& uniprot_to_interpro_table)
    {
        std::ifstream file{ uniprot_to_interpro_table };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " + uniprot_to_interpro_table.string() + ", aborting...");
        }

        std::string line;
        mipfinder::interpro::ProteinDomains results;
        while (std::getline(file, line)) {
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