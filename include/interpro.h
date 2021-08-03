#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>
 
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
        repeat,
        unknown
    };

    struct Entry
    {
        std::string description;
        Type type;
        auto operator<=>(const Entry&) const = default;
    };

    using DomainIdentifier = std::string;
    using Entries = std::unordered_map<DomainIdentifier, Entry>;
    using ProteinDomains = std::unordered_map<std::string, std::unordered_set<std::string>>;

    /* @interpro_entry_list is a tsv-file that specifies which type each InterPro identifier is
     * Column 1: InterPro identifier
     * Column 2: Identifier type
     * Column 3: Identifier description
     * 
     * Return an associative array where the keys are InterPro entry identifiers, and the values contain
     * data about the entry. See 'Data' struct for more information.
     */
    Entries parseEntryList(const std::filesystem::path& interpro_entry_list);

    //Map every UniProt ID to all the InterPro domains the protein contains.
    //
    //@uniprot_to_interpro_table is a tsv-file that specifies which UniProt entry has which InterPro
    //domains
    //Column 1: UniProt ID (without sequence version)
    //Column 2: Sequence version of that corresponding UniProt ID
    //Column 3: A comma (;) separated list of InterPro entry identifiers
    //
    //Return an associative array where the keys are UniProt identifiers (including the sequence version) and the
    //values are a set of InterPro identifiers corresponding to the InterPro entires that protein has been
    //annotated to have. If the input file is not in a correct format, the behaviour is unspecified.
    ProteinDomains parseProteinDomainList(const std::filesystem::path& uniprot_to_interpro_table);
}
#endif
