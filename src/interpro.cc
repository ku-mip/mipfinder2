#include "helpers.h"
#include "interpro.h"
#include "protein.h"

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
    std::vector<mipfinder::Interpro::Entry> parseInterproEntryList(const std::filesystem::path& interpro_entry_list)
    {
        std::ifstream file{ interpro_entry_list };
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open " + interpro_entry_list.string() +
                ", aborting...");
        }

        std::vector<mipfinder::Interpro::Entry> results;
        std::string line;
        std::getline(file, line); //Skip the header in database file
        while (std::getline(file, line)) {
            const auto tokens = mipfinder::tokenise(line, '\t');

            //Correctly formatted file has three columns of data
            if (tokens.size() != 3) {
                continue;
            }

            mipfinder::Interpro::Entry::Type type;
            const auto& domain_type = tokens[1];
            if (domain_type == "Active_site") {
                type = mipfinder::Interpro::Entry::Type::active_site;
            }
            else if (domain_type == "Binding_site") {
                type = mipfinder::Interpro::Entry::Type::binding_site;
            }
            else if (domain_type == "Conserved_site") {
                type = mipfinder::Interpro::Entry::Type::conserved_site;
            }
            else if (domain_type == "Domain") {
                type = mipfinder::Interpro::Entry::Type::domain_type;
            }
            else if (domain_type == "Family") {
                type = mipfinder::Interpro::Entry::Type::family;
            }
            else if (domain_type == "Homologous_superfamily") {
                type = mipfinder::Interpro::Entry::Type::homologous_superfamily;
            }
            else if (domain_type == "PTM") {
                type = mipfinder::Interpro::Entry::Type::ptm;
            }
            else if (domain_type == "Repeat") {
                type = mipfinder::Interpro::Entry::Type::repeat;
            }
            else {
                type = mipfinder::Interpro::Entry::Type::unknown;
            }

            const auto& interpro_accession = tokens[0];
            const auto& entry_name = tokens[2];
            results.emplace_back(mipfinder::Interpro::Entry{
                .name = entry_name, .accession = interpro_accession, .type = type });
        }
        return results;
    }
} // namespace detail

namespace mipfinder
{
    Interpro::Interpro(const std::filesystem::path& interpro_entry_list)
        : m_entries(detail::parseInterproEntryList(interpro_entry_list))
    {
        std::sort(std::begin(m_entries), std::end(m_entries));
    }

    Interpro::const_iterator Interpro::begin() const
    {
        return m_entries.begin();
    }

    Interpro::const_iterator Interpro::cbegin() const
    {
        return m_entries.cbegin();
    }

    Interpro::const_iterator Interpro::end() const
    {
        return m_entries.end();
    }

    Interpro::const_iterator Interpro::cend() const
    {
        return m_entries.cend();
    }

    Interpro::const_iterator find(const Interpro& interpro_database, const std::string& entry_accession)
    {
        auto comparator = [](const Interpro::Entry& element,
            const std::string& accession_to_find) {
                return element.accession < accession_to_find;
        };

        auto returned_entry = std::lower_bound(
            std::cbegin(interpro_database), std::cend(interpro_database), entry_accession, comparator);

        if (returned_entry != std::cend(interpro_database) &&
            (*returned_entry).name == entry_accession) {
            return returned_entry;
        }
        else {
            return interpro_database.cend();
        }
    }

} // namespace mipfinder