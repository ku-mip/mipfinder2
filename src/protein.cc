#include <array>
#include <cmath>
#include <regex>
#include <unordered_map>

#include "fasta.h"
#include "protein.h"

namespace detail
{
	/**
	 * @brief Removes all illegal characters from a given protein sequence
	 */
	std::string sanitiseProteinSequenceInput(std::string input_sequence)
	{
		static constexpr std::array<unsigned char, 21> allowed_amino_acid_characters {
			'A', 'C', 'D', 'E', 'F',
			'G', 'H', 'I', 'K', 'L',
			'M', 'N', 'P', 'Q', 'R',
			'S', 'T', 'V', 'W', 'Y',
			'*'
		};

		std::string sanitised_string;
		for (const char& original_character : input_sequence) {
			if (std::find(std::begin(allowed_amino_acid_characters),
						  std::end(allowed_amino_acid_characters),
						  std::toupper(original_character)) != std::end(allowed_amino_acid_characters)) {
				sanitised_string += original_character;
			}
		}
		return sanitised_string;
	}

	struct UniProtHeader
	{
		unsigned int sequence_version;
		unsigned int existence_level;
		std::string accession;
		std::string description;
	};

	/* Takes a UniProt database FASTA header as input and returns an array with the
	 * following elements:
	 * 1: UniProt accession name
	 * 2: Entry sequence version
	 * 3: Entry description
	 * 4: Entry protein existence level
	 *
	 * If the header cannot be successfully parsed, returns an empty array.
	 */
	UniProtHeader extractUniprotHeader(const std::string& header)
	{
		/* This regex currently only grabs the accession name, description and
		 * sequence version because the other data is not interesting */
		static const std::regex uniprot_regex(R"((?:>)(\w+)(?:\|)(\w+)(?:\|)(\w+)[ ](.+?)[ ](?:OS=).+(?:PE=)([0-9]).+(?:SV=)([0-9]))");

		std::string accession_name;
		std::string sequence_version;
		std::string description;
		std::string existence_level;

		std::smatch matches;
		if (std::regex_search(header, matches, uniprot_regex)) {

			if (matches.size() != 7) {
				return UniProtHeader{};
			}

			/* For a standard UniProt header the matches will as following:
			 * Match 0 - whole match
			 * Match 1 - database type
			 * Match 2 - UniProt accession name
			 * Match 3 - Entry name
			 * Match 4 - Description
			 * Match 5 - Protein existence level
			 * Match 6 - Sequence version */
			accession_name = matches[2];
			description = matches[4];
			existence_level = matches[5];
			sequence_version = matches[6];

			return UniProtHeader{ .sequence_version = std::stoul(sequence_version),
					 .existence_level = std::stoul(existence_level),
					 .accession = accession_name,
					 .description = description };
		}
		return UniProtHeader{};
	}
}


namespace mipfinder::protein
{	
	Identifier::Identifier(std::string identifier, unsigned int sequence_version)
		: m_identifier(identifier + Identifier::identifier_delimiter + std::to_string(sequence_version)) {}

	std::string Identifier::to_string() const
	{
		return m_identifier;
	}

	std::ostream& operator<<(std::ostream& os, const Identifier& obj)
	{
		os << obj.to_string();
		return os;
	}

	Sequence::Sequence(std::string sequence)
		: m_sequence(detail::sanitiseProteinSequenceInput(sequence)) { }

	std::size_t Sequence::length() const
	{
		return m_sequence.length();
	}

	std::string Sequence::to_string() const
	{
		return m_sequence;
	}

	Protein::Protein() : m_existence_level(0), m_type(Protein::Type::UNKNOWN)
	{
	}

	Protein::Protein(Identifier identifier,
					 Sequence sequence,
					 std::string description,
					 const unsigned int existence_level)
		: m_identifier(identifier),
		  m_sequence(sequence),
		  m_description(description),
		  m_existence_level(existence_level),
		  m_type(Protein::Type::UNKNOWN)
	{
	}


	Protein::Protein(const Protein& prot) :
		m_identifier(prot.m_identifier),
		m_sequence(prot.m_sequence),
		m_description(prot.m_description),
		m_existence_level(prot.m_existence_level),
		m_type(prot.m_type)
	{
	}

	std::ostream& operator<<(std::ostream& os, const Sequence& obj)
	{
		os << obj.to_string();
		return os;
	}

	void swap(Protein& lhs, Protein& rhs)
	{
		using std::swap;
		swap(lhs.m_identifier, rhs.m_identifier);
		swap(lhs.m_sequence, rhs.m_sequence);
		swap(lhs.m_description, rhs.m_description);
		swap(lhs.m_existence_level, rhs.m_existence_level);
		swap(lhs.m_type, rhs.m_type);
	}

	//Copy assignment operator
	Protein& Protein::operator=(Protein other)
	{
		swap(*this, other);
		return *this;
	}

	//Move constructor
	Protein::Protein(Protein&& other) noexcept : Protein()
	{
		swap(*this, other);
	}

	Identifier Protein::identifier() const
	{
		return m_identifier;
	}

	Sequence Protein::sequence() const
	{
		return m_sequence;
	}

	std::string Protein::description() const
	{
		return m_description;
	}

	unsigned int Protein::existenceLevel() const
	{
		return m_existence_level;
	}

	Proteome::Proteome(const std::filesystem::path& fasta_file)
	{
		
	}

	std::size_t Proteome::size() const noexcept
	{
		return proteins.size();
	}

	Proteome::iterator Proteome::find(const Identifier& identifier)
	{
		return proteins.find(identifier);
	}

	Proteome::const_iterator Proteome::find(const Identifier& identifier) const
	{
		return proteins.find(identifier);
	}

	/**
	 * @brief  Extract protein data from a FASTA file.
	 * @param  fasta_file  Location of the FASTA file.
	 *
	 * @throws  std::runtime_error  If @a fasta_file cannot be opened.
	 * @return  A collection of Proteins sorted by their identifier.
	 *
	 * The FASTA headers in the given @a fasta_file need to be in UniProt
	 * format. If they are not, the behaviour is unspecified.
	 */
	void Proteome::parseProteins(const std::filesystem::path& fasta_file)
	{
		const mipfinder::fasta::Entries fasta_entries = mipfinder::fasta::parse(fasta_file);
		for (const auto& entry : fasta_entries) {
			//TODO: Check if it is a uniprot header or not, if not then use the most
			//basic constructor for Proteins.
			const auto& header = detail::extractUniprotHeader(entry.header);
			const mipfinder::protein::Identifier identifier{
				header.accession,
				header.sequence_version };

			proteins.emplace(mipfinder::protein::Protein{
				identifier,
				entry.sequence,
				header.description,
				header.existence_level });
		}
	}

	double instability_index(const std::string& sequence)
	{
		static std::unordered_map<char, std::unordered_map<char, double>> dimw{
		  {'A', { {'A', 1.0 },   {'C', 44.94},  {'E', 1.0},   {'D', -7.49},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', -7.49},
				  {'K', 1.0},    {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 1.0},    {'P', 20.26},  {'S', 1.0},   {'R', 1.0},
				  {'T', 1.0},    {'W', 1.0},    {'V', 1.0},   {'Y', 1.0}    }},
		  {'C', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 20.26},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', 33.60},
				  {'K', 1.0},    {'M', 33.60},  {'L', 20.26}, {'N', 1.0},
				  {'Q', -6.54},  {'P', 20.26},  {'S', 1.0},   {'R', 1.0},
				  {'T', 33.60},  {'W', 24.68},  {'V', -6.54}, {'Y', 1.0}    }},
		  {'E', { {'A', 1.0},    {'C', 44.94},  {'E', 33.60}, {'D', 20.26},
				  {'G', 1.0},    {'F', 1.0},    {'I', 20.26}, {'H', -6.54},
				  {'K', 1.0},    {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 20.26},  {'P', 20.26},  {'S', 20.26}, {'R', 1.0},
				  {'T', 1.0},    {'W', -14.03}, {'V', 1.0},   {'Y', 1.0}    }},
		  {'D', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', 1.0},    {'F', -6.54},  {'I', 1.0},   {'H', 1.0},
				  {'K', -7.49},  {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 1.0},    {'P', 1.0},    {'S', 20.26}, {'R', -6.54},
				  {'T', -14.03}, {'W', 1.0},    {'V', 1.0},   {'Y', 1.0}    }},
		  {'G', { {'A', -7.49},  {'C', 1.0},    {'E', -6.54}, {'D', 1.0},
				  {'G', 13.34},  {'F', 1.0},    {'I', -7.49}, {'H', 1.0},
				  {'K', -7.49},  {'M', 1.0},    {'L', 1.0},   {'N', -7.49},
				  {'Q', 1.0},    {'P', 1.0},    {'S', 1.0},   {'R', 1.0},
				  {'T', -7.49},  {'W', 13.34},  {'V', 1.0},   {'Y', -7.49}  }},
		  {'F', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 13.34},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', 1.0},
				  {'K', -14.03}, {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 1.0},    {'P', 20.26},  {'S', 1.0},   {'R', 1.0},
				  {'T', 1.0},    {'W', 1.0},    {'V', 1.0},   {'Y', 33.601} }},
		  {'I', { {'A', 1.0},    {'C', 1.0},    {'E', 44.94}, {'D', 1.0},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', 13.34},
				  {'K', -7.49},  {'M', 1.0},    {'L', 20.26}, {'N', 1.0},
				  {'Q', 1.0},    {'P', -1.88},  {'S', 1.0},   {'R', 1.0},
				  {'T', 1.0},    {'W', 1.0},    {'V', -7.49}, {'Y', 1.0}    }},
		  {'H', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', -9.37},  {'F', -9.37},  {'I', 44.94}, {'H', 1.0},
				  {'K', 24.68},  {'M', 1.0},    {'L', 1.0},   {'N', 24.68},
				  {'Q', 1.0},    {'P', -1.88},  {'S', 1.0},   {'R', 1.0},
				  {'T', -6.54},  {'W', -1.88},  {'V', 1.0},   {'Y', 44.94}  }},
		  {'K', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', -7.49},  {'F', 1.0},    {'I', -7.49}, {'H', 1.0},
				  {'K', 1.0},    {'M', 33.60},  {'L', -7.49}, {'N', 1.0},
				  {'Q', 24.64},  {'P', -6.54},  {'S', 1.0},   {'R', 33.60},
				  {'T', 1.0},    {'W', 1.0},    {'V', -7.49}, {'Y', 1.0}    }},
		  {'M', { {'A', 13.34},  {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', 58.28},
				  {'K', 1.0},    {'M', -1.88},  {'L', 1.0},   {'N', 1.0},
				  {'Q', -6.54},  {'P', 44.94},  {'S', 44.94}, {'R', -6.54},
				  {'T', -1.88},  {'W', 1.0},    {'V', 1.0},   {'Y', 24.68}  }},
		  {'L', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', 1.0},
				  {'K', -7.49},  {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 33.60},  {'P', 20.26},  {'S', 1.0},   {'R', 20.26},
				  {'T', 1.0},    {'W', 24.68},  {'V', 1.0},   {'Y', 1.0}    }},
		  {'N', { {'A', 1.0},    {'C', -1.88},  {'E', 1.0},   {'D', 1.0},
				  {'G', -14.03}, {'F', -14.03}, {'I', 44.94}, {'H', 1.0},
				  {'K', 24.68},  {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', -6.54},  {'P', -1.88},  {'S', 1.0},   {'R', 1.0},
				  {'T', -7.49},  {'W', -9.37},  {'V', 1.0},   {'Y', 1.0}    }},
		  {'Q', { {'A', 1.0},    {'C', -6.54},  {'E', 20.26}, {'D', 20.26},
				  {'G', 1.0},    {'F', -6.54},  {'I', 1.0},   {'H', 1.0},
				  {'K', 1.0},    {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 20.26},  {'P', 20.26},  {'S', 44.94}, {'R', 1.0},
				  {'T', 1.0},    {'W', 1.0},    {'V', -6.54}, {'Y', -6.54}  }},
		  {'P', { {'A', 20.26},  {'C', -6.54},  {'E', 18.38}, {'D', -6.54},
				  {'G', 1.0},    {'F', 20.26},  {'I', 1.0},   {'H', 1.0},
				  {'K', 1.0},    {'M', -6.54},  {'L', 1.0},   {'N', 1.0},
				  {'Q', 20.26},  {'P', 20.26},  {'S', 20.26}, {'R', -6.54},
				  {'T', 1.0},    {'W', -1.88},  {'V', 20.26}, {'Y', 1.0}    }},
		  {'S', { {'A', 1.0},    {'C', 33.60},  {'E', 20.26}, {'D', 1.0},
				  {'G', 1.0},    {'F', 1.0},    {'I', 1.0},   {'H', 1.0},
				  {'K', 1.0},    {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 20.26},  {'P', 44.94},  {'S', 20.26}, {'R', 20.26},
				  {'T', 1.0},    {'W', 1.0},    {'V', 1.0},   {'Y', 1.0}    }},
		  {'R', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', -7.49},  {'F', 1.0},    {'I', 1.0},   {'H', 20.26},
				  {'K', 1.0},    {'M', 1.0},    {'L', 1.0},   {'N', 13.34},
				  {'Q', 20.26},  {'P', 20.26},  {'S', 44.94}, {'R', 58.28},
				  {'T', 1.0},    {'W', 58.28},  {'V', 1.0},   {'Y', -6.54}  }},
		  {'T', { {'A', 1.0},    {'C', 1.0},    {'E', 20.26}, {'D', 1.0},
				  {'G', -7.49},  {'F', 13.34},  {'I', 1.0},   {'H', 1.0},
				  {'K', 1.0},    {'M', 1.0},    {'L', 1.0},   {'N', -14.03},
				  {'Q', -6.54},  {'P', 1.0},    {'S', 1.0},   {'R', 1.0},
				  {'T', 1.0},    {'W', -14.03}, {'V', 1.0},   {'Y', 1.0}    }},
		  {'W', { {'A', -14.03}, {'C', 1.0},    {'E', 1.0},   {'D', 1.0},
				  {'G', -9.37},  {'F', 1.0},    {'I', 1.0},   {'H', 24.68},
				  {'K', 1.0},    {'M', 24.68},  {'L', 13.34}, {'N', 13.34},
				  {'Q', 1.0},    {'P', 1.0},    {'S', 1.0},   {'R', 1.0},
				  {'T', -14.03}, {'W', 1.0},    {'V', -7.49}, {'Y', 1.0}    }},
		  {'V', { {'A', 1.0},    {'C', 1.0},    {'E', 1.0},   {'D', -14.03},
				  {'G', -7.49},  {'F', 1.0},    {'I', 1.0},   {'H', 1.0},
				  {'K', -1.88},  {'M', 1.0},    {'L', 1.0},   {'N', 1.0},
				  {'Q', 1.0},    {'P', 20.26},  {'S', 1.0},   {'R', 1.0},
				  {'T', -7.49},  {'W', 1.0},    {'V', 1.0},   {'Y', -6.54}  }},
		  {'Y', { {'A', 24.68},  {'C', 1.0},    {'E', -6.54}, {'D', 24.68},
				  {'G', -7.49},  {'F', 1.0},    {'I', 1.0},   {'H', 13.34},
				  {'K', 1.0},    {'M', 44.94},  {'L', 1.0},   {'N', 1.0},
				  {'Q', 1.0},    {'P', 13.34},  {'S', 1.0},   {'R', -15.91},
				  {'T', -7.49},  {'W', -9.37},  {'V', 1.0},   {'Y', 13.34}  }}};

		double score{0};
		for (std::size_t pos = 0; pos != sequence.size() - 1; ++pos) {
			const auto dipeptide = sequence.substr(pos, 2);
			score += dimw[dipeptide[0]][dipeptide[1]];
		}
		return (10.0 / sequence.length()) * score;
	}

}