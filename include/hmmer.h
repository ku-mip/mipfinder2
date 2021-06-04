#ifndef MIPFINDER_HMMER_H
#define MIPFINDER_HMMER_H

#include <filesystem>
#include <functional>

#include "aliases.h"
#include "interpro.h"
#include "helpers.h"
#include "protein.h"

namespace mipfinder
{
	class Proteome;
	class Configuration;

	template <typename T, typename U>
	struct Homologue
	{
		T* query;
		U* target;
		double bitscore;
	};
}

namespace std
{
	template <typename T, typename U>
	struct hash<mipfinder::Homologue<T, U>>
	{
		std::size_t operator()(const mipfinder::Homologue<T, U>& k) const
		{
			using std::hash;

			return ((hash<T>()(k)
					^ (hash<U>()(k) << 1)) >> 1);
		}
	};
}

namespace mipfinder::hmmer
{
	struct Result;
	typedef std::vector<mipfinder::hmmer::Result> Results;

	//Runs phmmer with the query against the database. Extra phmmer options can be specified in the 
	//`extra_parameters` variable. 
	template <typename T, typename U>
	void phmmer(T& query,
				U& database,
				const std::filesystem::path& results_file,
				const std::string& extra_parameters)
	{
		//Convert query and database into FASTA files
		const std::filesystem::path results_path = results_file.parent_path();

		const std::filesystem::path query_fasta_file{"hmmer_query.fasta"};
		const std::filesystem::path query_file_location =
			results_path / query_fasta_file;

		const std::filesystem::path database_fasta_file{"hmmer_database.fasta"};
		const std::filesystem::path database_file_location =
			results_path / database_fasta_file;

		mipfinder::proteinToFasta(query, query_file_location);
		mipfinder::proteinToFasta(database, database_file_location);
		return phmmer(query_file_location, database_file_location, results_file, extra_parameters);
	}

	template <>
	void phmmer(const std::filesystem::path& query_file,
				const std::filesystem::path& database_file,
				const std::filesystem::path& results_file,
				const std::string& extra_parameters);


	void buildHmmerProfile(const std::filesystem::path& msa_file,
						   const std::filesystem::path& output_file,
						   const std::string& extra_parameters = "");

	void hmmsearch(const std::filesystem::path& profile_file,
				   const ProteinSet& database,
				   const std::filesystem::path& output_file,
				   const std::string& extra_parameters = "");

	void hmmsearch(const std::filesystem::path& profile_file,
				   const std::filesystem::path& database_file,
				   const std::filesystem::path& output_file,
				   const std::string& extra_parameters = "");

	/* Reads in a results file in table format (`--tblout` from HMMMER) */
	mipfinder::hmmer::Results parseResults(const std::filesystem::path& results_file);


	/* Keeps up to @hits_to_keep of best-scoring HMMER results for each query */
	mipfinder::hmmer::Results
		keepTopHits(const mipfinder::hmmer::Results& results, std::size_t hits_to_keep);

	/* Removes all entries from HMMER results files where the query is the same
	 * as the target */
	mipfinder::hmmer::Results
		removeSelfHits(const mipfinder::hmmer::Results& results);

	/* Returns a list of proteins from @proteins that match the identifiers found
	 * in HMMER @results */
	mipfinder::ProteinSet
		convertToProtein(const mipfinder::hmmer::Results& results,
						 const mipfinder::ProteinSet& proteins);

	//Creates HMMER profiles from all cMIP Multiple Sequence Alignments. 
	//Creates one HMMER profile per one MSA.
	void createHmmProfiles(const std::filesystem::path& msa_dir,
						   const std::filesystem::path& output_dir);

	//Creates one large HMMER profile file from all homologous cMIP HMMER
	//profiles. This simplifies downstream processing.
	std::filesystem::path
		concatenateHmmProfiles(const std::filesystem::path& hmmer_profile_dir,
							   const std::filesystem::path& output_file);

	/* Returns all HmmerResults queries where the target length - query length >=
	 * @min_length_difference */
	mipfinder::hmmer::Results
		filterByLengthDifference(const mipfinder::hmmer::Results& results,
								 const mipfinder::ProteinSet& proteins,
								 unsigned int min_length_difference);

	/* Takes HmmerResults and filters out all targets that do not contain the same
	 * protein family identifier as the query */
	mipfinder::hmmer::Results
		filterByProteinFamilyAndDomain(const mipfinder::hmmer::Results& results,
									   const mipfinder::Proteome& proteome);

	/* Returns all HmmerResults queries that have less or equal than
	 * @maximum_homologues targets */
	mipfinder::hmmer::Results
		filterByHomologueCount(const mipfinder::hmmer::Results& results,
							   unsigned int maximum_homologues);

	template <typename Comparator = std::greater_equal<> >
	mipfinder::hmmer::Results filterByBitscore(const mipfinder::hmmer::Results& input,
											 double bitscore_cutoff,
											 Comparator comp = Comparator())
	{
		mipfinder::hmmer::Results tmp;
		for (const auto& result : input) {
			if (result.query == result.target) {
				tmp.push_back(result);
			}

			else if (comp(result.bitscore, bitscore_cutoff)) {
				tmp.push_back(result);
			}
		}
		return tmp;
	}

}
#endif
