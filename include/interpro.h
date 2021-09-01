#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <filesystem>
#include <string>
#include <vector>

namespace mipfinder
{
    /**
     *  @class  Interpro interpro.h "include/interpro.h"
     *  @brief  Provide a collection of InterPro entries.
     *
     *  @invariant  1. The collection is always sorted according to the InterPro
     *              entry accession value in a lexicographical order.
     *              2. The entry accessions are unique.
     */
    class Interpro {
    public:
        /**
         *  @class  Entry interpro.h "include/interpro.h"
         *  @brief  Represents a single InterPro entry.
         *
         *  Every InterPro entry has an accession identifier, a name (also
         *  known as description), and an entry type.
         */
        struct Entry {
            enum class Type {
                active_site,
                binding_site,
                conserved_site,
                domain_type,
                family,
                homologous_superfamily,
                ptm,
                repeat,
                unknown
            };
            std::string name;
            std::string accession;
            Type type;
            auto operator<=>(const Entry&) const = default;
        };

        /**
         *  @brief  Initialise a collection of InterPro entries.
         *  @param  entries  A three-column tsv-file that details each InterPro
         *                   entry (available from
                             https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list).
         *                   Column 1: InterPro entry accession.
         *                   Column 2: InterPro entry type.
         *                   Column 3: InterPro entry name (description).
         *  @throw  std::runtime_error If @a entries cannot be opened.
         */
        Interpro(const std::filesystem::path& entries);

        using const_iterator = std::vector<Entry>::const_iterator;
        const_iterator begin() const;
        const_iterator cbegin() const;
        const_iterator end() const;
        const_iterator cend() const;

    private:
        std::vector<Entry> m_entries;
    };

    /**
     * @brief  Find InterPro entry by its accession.
     * @param  interpro_database  Database to search.
     * @param  entry_accession  Accession of the entry to search for.
     * @return  Read-only (constant) iterator to the given InterPro entry, or
     *          iterator equivalent to Interpro::cend() if the entry could not be
     *          be found.
     */
    Interpro::const_iterator find(const Interpro& interpro_database, const std::string& entry_accession);


    //using ProteinDomains = std::unordered_map<std::string, std::unordered_set<std::string>>;

    ////Map every UniProt ID to all the InterPro domains the protein contains.
    ////
    ////@uniprot_to_interpro_table is a tsv-file that specifies which UniProt entry has which InterPro
    ////domains
    ////Column 1: UniProt ID (without sequence version)
    ////Column 2: Sequence version of that corresponding UniProt ID
    ////Column 3: A comma (;) separated list of InterPro entry identifiers
    ////
    ////Return an associative array where the keys are UniProt identifiers (including the sequence version) and the
    ////values are a set of InterPro identifiers corresponding to the InterPro entires that protein has been
    ////annotated to have. If the input file is not in a correct format, the behaviour is unspecified.
    //ProteinDomains parseProteinDomainList(const std::filesystem::path& uniprot_to_interpro_table);
}
#endif
