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
	template <typename T>
	requires std::convertible_to<T, std::filesystem::path>
		mipfinder::HmmerResults parseResultsFile(const T& results_file)
	{
		auto stream = mipfinder::file::open(results_file);

		mipfinder::HmmerResults results;
		std::string line;
		while (std::getline(stream, line)) {
			if (line.front() == '#') { //'#' lines are comments
				continue;
			}

			const auto tokens = mipfinder::tokenise(line, ' ');
			const std::string query = tokens[2];
			const std::string target = tokens[0];
			double bitscore = stod(tokens[5]);

			results.push_back(mipfinder::hmmer::Result{query, target, bitscore});
		}
		return results;
	}

	/* Keeps up to @hits_to_keep of best-scoring HMMER results for each query */
	mipfinder::HmmerResults
		keepTopHits(const mipfinder::HmmerResults& results, std::size_t hits_to_keep);

	/* Removes all entries from HMMER results files where the query is the same
	 * as the target */
	mipfinder::HmmerResults
		removeSelfHits(const mipfinder::HmmerResults& results);

	/* Returns a list of proteins from @proteins that match the identifiers found
	 * in HMMER @results */
	mipfinder::ProteinSet
		convertToProtein(const mipfinder::HmmerResults& results,
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
	mipfinder::HmmerResults
		filterByLengthDifference(const mipfinder::HmmerResults& results,
								 const mipfinder::ProteinSet& proteins,
								 unsigned int min_length_difference);

	/* Takes HmmerResults and filters out all targets that do not contain the same
	 * protein family identifier as the query */
	mipfinder::HmmerResults
		filterByProteinFamilyAndDomain(const mipfinder::HmmerResults& results,
									   const mipfinder::Proteome& proteome);

	/* Returns all HmmerResults queries that have less or equal than
	 * @maximum_homologues targets */
	mipfinder::HmmerResults
		filterByHomologueCount(const mipfinder::HmmerResults& results,
							   unsigned int maximum_homologues);

	template <typename Comparator = std::greater_equal<> >
	mipfinder::HmmerResults filterByBitscore(const mipfinder::HmmerResults& input,
											 double bitscore_cutoff,
											 Comparator comp = Comparator())
	{
		mipfinder::HmmerResults tmp;
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
