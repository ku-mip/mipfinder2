#ifndef MIPFINDER_PROTEIN_H
#define MIPFINDER_PROTEIN_H

#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>

#include "ancestor.h"
#include "aliases.h"
#include "homologue.h"
#include "go.h"
#include "interpro.h"

namespace mipfinder
{
	using ProteinList = std::vector<Protein>;

	//Protein holds information about one specific protein isoform. The Protein identifier is composed
	//of UniProt accession name plus the sequence version.
	class Protein
	{
	public:
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
			m_existence_level(prot.m_existence_level), m_score(prot.m_score), m_type(prot.m_type)/*,
			ancestors_(prot.ancestors_), interpro_entries_(prot.interpro_entries_), go_entries_(prot.go_entries_)*/ { }

		friend void swap(Protein& lhs, Protein& rhs);
		Protein& operator=(Protein other);
		Protein(Protein&& other) noexcept;

		auto operator<=>(const Protein&) const = default;

		std::string identifier() const;
		std::string sequence() const;
		std::string description() const;

		//std::vector<mipfinder::Ancestor> ancestors() const;
		//void addAncestor(const Ancestor& ancestor);

		//std::vector<mipfinder::Result> homologues() const;
		//void addHomologue(const Result& homologue);

		Type type() const;
		//std::string type_to_string() const;
		//void setType(const mipfinder::Protein::Type& type);

		//mipfinder::Go::Entries goEntries() const;
		//void addGoEntry(const mipfinder::Go::Entry& entry);

		//mipfinder::Interpro::Entries interproEntries() const;
		//void addInterproEntry(const mipfinder::Interpro::Entry& entry);

		int existenceLevel() const;

		std::size_t length() const;
		double score() const;
		void changeScore(double score);
	private:
		Protein() = default;

		std::string m_identifier;
		std::string m_sequence;
		std::string m_description;
		int m_existence_level;
		double m_score;

		mipfinder::Protein::Type m_type;
		//std::vector<mipfinder::Ancestor> ancestors_;
		//std::vector<mipfinder::Result> homologues_;

		//mipfinder::Interpro::Entries interpro_entries_;
		//mipfinder::Go::Entries go_entries_;
	};

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