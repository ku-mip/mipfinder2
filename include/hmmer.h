#ifndef MIPFINDER_HMMER_H
#define MIPFINDER_HMMER_H

#include <filesystem>
#include <functional>

#include "aliases.h"
#include "interpro.h"
#include "helpers.h"

namespace mipfinder
{
	class Proteome;
	class Configuration;
}

namespace mipfinder::hmmer
{
	struct Result;
	typedef std::vector<mipfinder::hmmer::Result> Results;

	//Runs phmmer with the query against the database. Extra phmmer options can be specified in the 
	//`extra_parameters` variable. 
	void phmmer(const ProteinSet& query,
				const ProteinSet& database,
				const std::filesystem::path& output_file,
				const std::string& extra_parameters = "");

	void phmmer(const std::filesystem::path& query_file,
				const std::filesystem::path& database_file,
				const std::filesystem::path& output_file,
				const std::string& extra_parameters = "");

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
