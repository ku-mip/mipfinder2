#ifndef MIPFINDER_INTERPRO_H
#define MIPFINDER_INTERPRO_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>
 
namespace mipfinder::interpro
{
    //Entry class denotes a unique InterPro entry. Every possible InterPro annotation for a protein
    //is composed from InterPro entries. 
    struct Entry
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
        std::string name;
        std::string accession;
        Type type;
        auto operator<=>(const Entry&) const = default;
    };

    using Entries = std::vector<Entry>;
    /* @interpro_entry_list is a tsv-file that details each InterPro entry (available from https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list [accessed 06/08/21]).
     * Column 1: InterPro identifier
     * Column 2: Identifier type
     * Column 3: Identifier description
     * 
     * Exceptions
     *  If `interpro_entry_list` cannot be opened, throws std::runtime_error
     * 
     * Return a list of InterPro entries in indeterminate order.
     */
    Entries parseEntryList(const std::filesystem::path& interpro_entry_list);

    using ProteinDomains = std::unordered_map<std::string, std::unordered_set<std::string>>;

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
