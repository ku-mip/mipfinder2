#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <filesystem>
#include <string>
#include <vector>
#include <unordered_map>

#include "protein.h"

namespace mipfinder::interpro
{
    class IdentifierMapping
    {
    public:
        using InterproIdentifier = std::string;
        using Domains = std::vector<InterproIdentifier>;
        using DomainTable = std::unordered_map<mipfinder::protein::Identifier, Domains>;

        Domains domains(const mipfinder::protein::Identifier& identifier);

        IdentifierMapping(const std::filesystem::path& identifier_to_domains);
    private:

        DomainTable identifier_domains;
    };


    /**
     *  @class  Interpro interpro.h "include/interpro.h"
     *  @brief  Provide a collection of InterPro entries.
     *
     *  @invariant  1. The collection is always sorted according to the InterPro
     *              entry accession value in a lexicographical order.
     *              2. The entry accessions are unique.
     */
    class Database {
    public:
        /**
         *  @class  Entry interpro.h "include/interpro.h"
         *  @brief  Represents a single InterPro entry.
         *
         *  Every InterPro entry has an accession identifier, a name (also
         *  known as description), and an entry type.
         */
        struct Entry {
            /**
             * @brief  Denotes the type of an InterPro entry annotation.
             */
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
        Database(const std::filesystem::path& entry_file);

        using Entries = std::vector<Entry>;
        using const_iterator = Entries::const_iterator;

        const_iterator begin() const;
        const_iterator cbegin() const;

        const_iterator end() const;
        const_iterator cend() const;
    private:
        Entries entries;
    };

    /**
     * @brief  Find InterPro entry by its accession.
     * @param  interpro_database  Database to search.
     * @param  entry_accession  Accession of the entry to search for.
     * @return  Read-only (constant) iterator to the given InterPro entry, or
     *          iterator equivalent to Interpro::cend() if the entry could not be
     *          be found.
     */
    Database::const_iterator find(const Database& interpro_database, const std::string& entry_accession);
}
#endif
