#include <algorithm>

#include "helpers.h"
#include "interpro.h"

namespace detail {

    /* Return a list of individual InterPro entries.
     *
     * @interpro_entry_list is a tsv-file that details each InterPro entry (available from https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list [accessed 06/08/21]).
     * Column 1: InterPro identifier
     * Column 2: Identifier type
     * Column 3: Identifier description
     *
     * The first line of the file must contain the headers for the columns. If the
     * file does not correspond to the given format, the behaviour is undefined.
     *
     * Exceptions
     *  If `interpro_entry_list` cannot be opened, throws std::runtime_error
     *
     * Return a list of InterPro entries in indeterminate order.
     */
    mipfinder::interpro::Database::Entries parseInterproEntryList(const std::filesystem::path& interpro_entry_list)
    {
        std::ifstream file{ interpro_entry_list };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " + interpro_entry_list.string() +
                ", aborting...");
        }

        mipfinder::interpro::Database::Entries results;
        std::string line;
        std::getline(file, line); //Skip the header in database file
        while (std::getline(file, line)) {
            const auto tokens = mipfinder::tokenise(line, '\t');

            //Correctly formatted file has three columns of data
            if (tokens.size() != 3) {
                continue;
            }

            mipfinder::interpro::Database::Entry::Type type;
            const auto& domain_type = tokens[1];
            if (domain_type == "Active_site") {
                type = mipfinder::interpro::Database::Entry::Type::active_site;
            }
            else if (domain_type == "Binding_site") {
                type = mipfinder::interpro::Database::Entry::Type::binding_site;
            }
            else if (domain_type == "Conserved_site") {
                type = mipfinder::interpro::Database::Entry::Type::conserved_site;
            }
            else if (domain_type == "Domain") {
                type = mipfinder::interpro::Database::Entry::Type::domain_type;
            }
            else if (domain_type == "Family") {
                type = mipfinder::interpro::Database::Entry::Type::family;
            }
            else if (domain_type == "Homologous_superfamily") {
                type = mipfinder::interpro::Database::Entry::Type::homologous_superfamily;
            }
            else if (domain_type == "PTM") {
                type = mipfinder::interpro::Database::Entry::Type::ptm;
            }
            else if (domain_type == "Repeat") {
                type = mipfinder::interpro::Database::Entry::Type::repeat;
            }
            else {
                type = mipfinder::interpro::Database::Entry::Type::unknown;
            }

            const auto& interpro_accession = tokens[0];
            const auto& entry_name = tokens[2];
            results.emplace_back(mipfinder::interpro::Database::Entry{
                .name = entry_name, .accession = interpro_accession, .type = type });
        }
        return results;
    }

    /**
     * @brief  Associate each UniProt identifier with their InterPro domains
     * @param  mapping_file  A tsv-file that specifies which UniProt entry has which InterPro
     *                       domains.
     *                       Column 1: UniProt accession without sequence version.
     *                       Column 2: Sequence version of the corresponding UniProt accession.
     *                       Column 3: Commma (;) separated list of InterPro entry identifiers.
     * @return  Associative array where each key is a protein identifier and
     *          the values are a list of unique InterPro identifiers.
     *
     * If the input file is not in a correct format, the behaviour is unspecified.
     */
    mipfinder::interpro::IdentifierMapping::DomainTable
    parseProteinDomainList(const std::filesystem::path& mapping_file)
    {
        std::ifstream file{ mapping_file };
        if (!file.is_open()) {
            return mipfinder::interpro::IdentifierMapping::DomainTable{};
        }
        mipfinder::interpro::IdentifierMapping::DomainTable parsed_list{};
        std::string line;
        while (std::getline(file, line)) {
            auto tokens = mipfinder::tokenise(line, '\t');
            std::string uniprot_accession = tokens[0];
            unsigned int sequence_version = std::stoul(tokens[1]);
            mipfinder::protein::Identifier id{ uniprot_accession, sequence_version };
            auto interpro_identifiers = mipfinder::tokenise(tokens[2], ';');

            for (const auto& interpro_id : interpro_identifiers) {
                parsed_list[id].push_back(interpro_id);
            }
        }
        return parsed_list;
    }

} // namespace detail

namespace mipfinder::interpro
{
    IdentifierMapping::IdentifierMapping(const std::filesystem::path& identifier_to_domains) : identifier_domains(detail::parseProteinDomainList(identifier_to_domains)) { }

    IdentifierMapping::Domains IdentifierMapping::domains(const mipfinder::protein::Identifier& identifier)
    {
        return identifier_domains.contains(identifier) ? identifier_domains.at(identifier) : IdentifierMapping::Domains{};
    }


    Database::Database(const std::filesystem::path& entries)
        : entries(detail::parseInterproEntryList(entries))
    {
        std::sort(std::begin(entries), std::end(entries));
    }

    Database::const_iterator Database::begin() const
    {
        return entries.begin();
    }

    Database::const_iterator Database::cbegin() const
    {
        return entries.cbegin();
    }

    Database::const_iterator Database::end() const
    {
        return entries.end();
    }

    Database::const_iterator Database::cend() const
    {
        return entries.cend();
    }

    Database::const_iterator find(const Database& interpro_database, const std::string& entry_accession)
    {
        auto comparator = [](const auto& element,
            const std::string& accession_to_find) {
                return element.accession < accession_to_find;
        };

        auto returned_entry = std::lower_bound(
            std::cbegin(interpro_database), std::cend(interpro_database), entry_accession, comparator);

        if (returned_entry != std::cend(interpro_database) &&
            (*returned_entry).name == entry_accession) {
            return returned_entry;
        }

        return interpro_database.cend();
    }

} // namespace mipfinder