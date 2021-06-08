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
	//Constraints
	template <typename T>
	concept ProteinRange = std::ranges::range<T> && std::same_as<typename std::ranges::range_value_t<T>, mipfinder::Protein>;

	template <typename T>
	concept HasIdentifier = requires (T t)
	{
		{ t.identifier() } -> std::convertible_to<std::string>;
	};

	template <typename T>
	concept HasSequence = requires (T t)
	{
		{ t.sequence() } -> std::convertible_to<std::string>;
	};

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


		//Constructors

		Protein(const std::string& identifier,
				const std::string& sequence,
				const std::string& description,
				int protein_existence);
		~Protein() = default;

		//Copy constructor
		Protein(const Protein& prot) :
			identifier_(prot.identifier_), sequence_(prot.sequence_), description_(prot.description_),
			existence_level_(prot.existence_level_), score_(prot.score_), type_(prot.type_),
			ancestors_(prot.ancestors_), interpro_entries_(prot.interpro_entries_), go_entries_(prot.go_entries_) { }

		//Copy assignment operator
		Protein& operator=(Protein other)
		{
			std::swap(*this, other);

			//identifier_ = std::move(other.identifier_);
			//sequence_ = std::move(other.sequence_);
			//description_ = std::move(other.description_);
			//existence_level_ = std::move(other.existence_level_);
			//score_ = std::move(other.score_);

			//type_ = std::move(other.type_);
			//ancestors_ = std::move(other.ancestors_);
			//interpro_entries_ = std::move(other.interpro_entries_);
			//go_entries_ = std::move(go_entries_);
			return *this;
		}


		//Move constructor
		Protein(Protein&& other) : Protein()
		{
			std::swap(*this, other);
		}

		auto operator<=>(const Protein& protein) const = default;

		std::string identifier() const;
		std::string sequence() const;
		std::string description() const;

		std::vector<mipfinder::Ancestor> ancestors() const;
		void addAncestor(const Ancestor& ancestor);

		//std::vector<mipfinder::Result> homologues() const;
		//void addHomologue(const Result& homologue);

		Type type() const;
		std::string type_to_string() const;
		void setType(const mipfinder::Protein::Type& type);

		mipfinder::Go::Entries goEntries() const;
		void addGoEntry(const mipfinder::Go::Entry& entry);

		mipfinder::Interpro::Entries interproEntries() const;
		void addInterproEntry(const mipfinder::Interpro::Entry& entry);

		int existenceLevel() const;

		std::size_t length() const;
		double score() const;
		void changeScore(double score);
	private:
		Protein() = default;

		std::string identifier_;
		std::string sequence_;
		std::string description_;
		int existence_level_;
		double score_;

		mipfinder::Protein::Type type_;
		std::vector<mipfinder::Ancestor> ancestors_;
		//std::vector<mipfinder::Result> homologues_;

		mipfinder::Interpro::Entries interpro_entries_;
		mipfinder::Go::Entries go_entries_;
	};

	/* Returns the instability index of the sequence */
	double instability_index(const std::string& sequence);



	//bool operator==(const Protein& lhs, const Protein& rhs);

	//Creates a FASTA file from existing Proteins
	template <typename T>
	requires std::ranges::range<T>
		&& mipfinder::HasIdentifier<std::ranges::range_value_t<T>>
		&& mipfinder::HasSequence<std::ranges::range_value_t<T>>
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