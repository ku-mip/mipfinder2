#include <array>

#include <cmath>
#include <unordered_map>

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
}


namespace mipfinder::protein
{
	Identifier::Identifier() : m_identifier(std::string{}) {}
	
	Identifier::Identifier(std::string identifier, unsigned int sequence_version)
		: m_identifier(identifier +
					   Identifier::identifier_delimiter +
					   std::to_string(sequence_version))
	{
	}


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


	//std::vector<Ancestor> Protein::ancestors() const
	//{
	//	return ancestors_;
	//}

	//void Protein::addAncestor(const Ancestor& ancestor)
	//{
	//	ancestors_.push_back(ancestor);
	//}

	//std::vector<Result> Protein::homologues() const
	//{
	//  return homologues_;
	//}

	//void Protein::addHomologue(const Result& homologue)
	//{
	//  homologues_.push_back(homologue);
	//}

	//Protein::Type Protein::type() const
	//{
	//	return m_type;
	//}


	//std::string Protein::type_to_string() const
	//{
	//	if (m_type == Protein::Type::SINGLE_COPY) {
	//		return "Single-copy";
	//	}
	//	else if (m_type == Protein::Type::HOMOLOGOUS) {
	//		return "Homologous";
	//	}
	//	else if (m_type == Protein::Type::KNOWN_MIP) {
	//		return "Known MIP";
	//	}
	//	else {
	//		return "UNKNOWN TYPE";
	//	}
	//}

	//void Protein::setType(const Protein::Type& type)
	//{
	//	m_type = type;
	//}

	//std::vector<Go::Entry> Protein::goEntries() const
	//{
	//	return go_entries_;
	//}

	//void Protein::addGoEntry(const Go::Entry& entry)
	//{
	//	go_entries_.push_back(entry);
	//}

	//Interpro::Entries Protein::interproEntries() const
	//{
	//	return interpro_entries_;
	//}

	//void Protein::addInterproEntry(const Interpro::Entry& entry)
	//{
	//	interpro_entries_.push_back(entry);
	//}



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