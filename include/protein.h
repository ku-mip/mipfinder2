#ifndef MIPFINDER_PROTEIN_H
#define MIPFINDER_PROTEIN_H

#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "hmmer.h"
#include "easylogging++.h"



namespace mipfinder::protein
{
    /**
     * @brief  Represents a protein identifier.
     */
    class Identifier
    {
    public:
        Identifier() = default;
        Identifier(std::string identifier, unsigned int sequence_version);
        std::string to_string() const;
        auto operator<=>(const Identifier& other) const = default;
    private:
        /* Used to delimit the protein identifier from its sequence version */
        static constexpr char identifier_delimiter = '-';
        std::string m_identifier;
    };
}

namespace std
{
    template <>
    struct hash<mipfinder::protein::Identifier>
    {
        std::size_t operator()(const mipfinder::protein::Identifier& k) const
        {
            using std::hash;
            return hash<std::string>()(k.to_string());
        }
    };
}

namespace mipfinder::protein
{
    /**
     * @brief  Writes the Identifier into an output stream.
     */
    std::ostream& operator<<(std::ostream& os, const Identifier& obj);

    /**
     * @brief  Represents a protein sequence in single-amino acid code.
     *
     * Uses the canonical NCBI single-amino acid lettering to represent all
     * canonical 20 amino acids.
     * (Available from https://www.ncbi.nlm.nih.gov/books/NBK20381/).
     * Stop codons are represented by asterisks (*).
     */
    class Sequence
    {
    public:
        Sequence() = default;

        /**
         * @brief  Constructs a sequence from the given string of amino acids.
         * @param  sequence  An amino acid sequence in single-character format.
         *
         * The constructor will sanitise the input and remove all disallowed
         * characters from the @a sequence.
         */
        Sequence(std::string sequence);
        /**
         * @brief  Converts the sequence into a string.
         */
        std::string to_string() const;
        /**
         * @brief  Returns the number of amino acids in the sequence.
         */
        std::size_t length() const;
    private:
        std::string m_sequence;
    };

    /**
     * @brief  Writes the Sequence into an output stream.
     */
    std::ostream& operator<<(std::ostream& os, const Sequence& obj);
}

namespace std
{
    template <>
    struct hash<mipfinder::protein::Sequence>
    {
        std::size_t operator()(const mipfinder::protein::Sequence& k) const
        {
            using std::hash;
            return hash<std::string>()(k.to_string());
        }
    };
}

namespace mipfinder::protein
{
    /**
     * @brief  Managing a protein
     */
    class Protein
    {
    public:
        /**
         * @brief  Declares the type of a protein. Possibly deprecated.
         */
        enum class Type
        {
            UNKNOWN,
            CMIP,
            SINGLE_COPY,
            HOMOLOGOUS,
            ANCESTOR,
            KNOWN_MIP
        };

        Protein(Identifier identifier,
            Sequence sequence,
            std::string description,
            unsigned int protein_existence);
        ~Protein() = default;
        Protein(const Protein& prot);
        Protein& operator=(Protein other);
        Protein(Protein&& other) noexcept;

        friend void swap(Protein& lhs, Protein& rhs);
        auto operator<=>(const Protein&) const = default;

        Identifier identifier() const;
        Sequence sequence() const;
        std::string description() const;
        unsigned int existenceLevel() const;


    private:
        Protein();

        Identifier m_identifier;
        Sequence m_sequence;
        std::string m_description;
        unsigned int m_existence_level;
        mipfinder::protein::Protein::Type m_type;
    };
}

namespace std
{
    template <>
    struct hash<mipfinder::protein::Protein>
    {
        std::size_t operator()(const mipfinder::protein::Protein& k) const
        {
            using std::hash;
            return ((hash<mipfinder::protein::Identifier>()(k.identifier())
                ^ (hash<mipfinder::protein::Sequence>()(k.sequence()) << 1)) >> 1);
        }
    };
}

namespace mipfinder::protein
{
    class Proteome
    {
    public:
        Proteome() = default;
        Proteome(const std::filesystem::path& fasta_file);
        
        using Proteins = std::unordered_map<Identifier, Protein>;
        using iterator = Proteins::iterator;
        using const_iterator = Proteins::const_iterator;
        
        iterator begin();
        const_iterator begin() const;
        const_iterator cbegin() const;

        iterator end();
        const_iterator end() const;
        const_iterator cend() const;

        iterator find(const Identifier& identifier);
        const_iterator find(const Identifier& identifier) const;

        std::size_t size() const noexcept;
    private:
        Proteins proteins;
        void parseProteins(const std::filesystem::path& fasta_file);
    };

    using ProteinList = std::vector<Protein>;

    /**
     * @brief  Calculates the instability index of the given protein sequence
     * @param  
     */
    double instability_index(const mipfinder::protein::Sequence& sequence);

    //template <typename T>
    //concept OstreamWriteable = requires(std::ofstream & out, T t)
    //{
    //    out << t;
    //};

    //TODO: Fix up these functions
    //Creates a FASTA file from a collection of proteins
    template <typename T>
    //requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> t)
    //{
    //    //{ t.identifier() } -> std::convertible_to<std::string>;
    //    //{ t.sequence() } -> std::convertible_to<std::string>;
    //}
    void proteinToFasta(T proteins, const std::filesystem::path& output)
    {
        std::ofstream f{ output };
        if (!f.is_open()) {
            return;
        }

        for (const auto& protein : proteins) {
            f << ">" << protein.identifier() << "\n" << protein.sequence() << "\n";
        }
    }

    template <typename T>
    void homologyToProtein(const T& proteome, const mipfinder::homology::Results& homology_search_results)
    {
        for (const auto& result : homology_search_results) {
            continue;
        }
    }

}
#endif
