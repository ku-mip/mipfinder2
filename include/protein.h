#ifndef MIPFINDER_PROTEIN_H
#define MIPFINDER_PROTEIN_H

#include <fstream>
#include <string>
#include <unordered_set>

#include "easylogging++.h"

namespace mipfinder
{
	//Protein holds information about one specific protein isoform. The Protein identifier is composed
	//of UniProt accession name plus the sequence version.
	class Protein
	{
	public:
		//This is used to delimit the protein UniProt identifier from its sequence version, e.g.
		//the final protein identifier is 'uniprot_id' + delimiter + 'sequence_version'.
		static constexpr char id_delimiter = '-'; 

		using Identifier = std::string;

		enum class Type
		{
			UNKNOWN,
			CMIP,
			SINGLE_COPY,
			HOMOLOGOUS,
			ANCESTOR,
			KNOWN_MIP
		};

		Protein(const std::string& identifier,
				const std::string& sequence,
				const std::string& description,
				int protein_existence);
		~Protein() = default;

		Protein(const Protein& prot) :
			m_identifier(prot.m_identifier), m_sequence(prot.m_sequence), m_description(prot.m_description),
			m_existence_level(prot.m_existence_level), m_type(prot.m_type) { }

		friend void swap(Protein& lhs, Protein& rhs);
		Protein& operator=(Protein other);
		Protein(Protein&& other) noexcept;

		auto operator<=>(const Protein&) const = default;

		std::string identifier() const;
		std::string sequence() const;
		std::string description() const;
		int existenceLevel() const;
		std::size_t length() const;
	private:
		Protein() = default;

		std::string m_identifier;
		std::string m_sequence;
		std::string m_description;
		int m_existence_level;
		mipfinder::Protein::Type m_type;
	};

	using ProteinList = std::vector<Protein>;

	/* Returns the instability index of the sequence */
	double instability_index(const std::string& sequence);

	//Creates a FASTA file from a collection of proteins
	template <typename T>
	requires std::ranges::range<T>&& requires (std::ranges::range_value_t<T> t)
	{
			{ t.identifier() } -> std::convertible_to<std::string>;
			{ t.sequence() } -> std::convertible_to<std::string>;
	}
	void proteinToFasta(T proteins, const std::filesystem::path& output)
	{
		LOG(DEBUG) << "Writing " << output << " FASTA file";
		std::ofstream f;
		f.open(output);
		for (const auto& protein : proteins) {
			f << ">" << protein.identifier() << "\n" << protein.sequence() << "\n";
		}
	}
}

namespace std
{
	template <>
	struct hash<mipfinder::Protein>
	{
		std::size_t operator()(const mipfinder::Protein& k) const
		{
			using std::hash;

			return ((hash<std::string>()(k.identifier())
					^ (hash<std::string>()(k.sequence()) << 1)) >> 1);
		}
	};
}


#endif