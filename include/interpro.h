#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <optional>
 
namespace mipfinder::interpro
{
    enum class Type {
        active_site,
        binding_site,
        conserved_site,
        domain_type,
        family,
        homologous_superfamily,
        ptm,
        repeat
    };

    struct Data
    {
        std::string description;
        Type type;
        auto operator<=>(const Data&) const = default;
    };

    using DomainIdentifier = std::string;
    using Entries = std::unordered_map<DomainIdentifier, Data>;
    using ProteinDomains = std::unordered_map < std::string, std::unordered_set<std::string>>;
    /* @database_file is a tsv-file that specifies which type each InterPro identifier is
     * Column 1: InterPro identifier
     * Column 2: Identifier type
     * Column 3: Identifier description
     */
    Entries parseEntryList(const std::filesystem::path& interpro_entry_list);
    ProteinDomains parseProteinDomainList(const std::filesystem::path& uniprot_to_interpro_table);
}
#endif
