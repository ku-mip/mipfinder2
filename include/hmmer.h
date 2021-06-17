#ifndef MIPFINDER_HMMER_H
#define MIPFINDER_HMMER_H

#include <filesystem>
#include <functional>
#include <limits>
#include <ranges>

#include "aliases.h"
#include "interpro.h"
#include "helpers.h"
#include "protein.h"

namespace mipfinder::homology
{
	struct Result
	{
		std::string query;
		std::string target;
		double bitscore;
	};

	using Results = std::vector<mipfinder::homology::Result>;

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
		const std::filesystem::path query_file_location = results_path / "hmmer_query.fasta";
		const std::filesystem::path database_file_location = results_path / "hmmer_database.fasta";

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
	mipfinder::homology::Results parseResults(const std::filesystem::path& results_file);

	//Find the corresponding proteins from the homology search results
	template <typename Cont>
	requires std::ranges::range<Cont>
	Cont findCorrespondingProteins(const mipfinder::homology::Results& results,
														 const Cont& proteome)
	{
		/* Lookup table for fast searching */
		std::unordered_map<std::string, mipfinder::Protein> lookup_table;
		for (const auto& protein : proteome) {
			lookup_table.insert(std::make_pair(protein.identifier(), protein));
		}

		Cont found_proteins;
		for (const auto& result : results) {
			if (lookup_table.contains(result.query)) {
				found_proteins.push_back(lookup_table.at(result.query));
			}
		}

		//Remove duplicates
		std::sort(std::begin(found_proteins), std::end(found_proteins));
		auto new_last_element = std::unique(std::begin(found_proteins), std::end(found_proteins));
		found_proteins.erase(new_last_element, found_proteins.end());
		return found_proteins;
	}


	template <typename T>
	requires std::ranges::range<T> && requires (std::ranges::range_value_t<T> t)
	{
		{ t.bitscore } -> std::convertible_to<double>;
	}
	auto filterByBitscore(const T& homology_results,
						  const double minimum_bitscore = (std::numeric_limits<double>::min)(), //Min and max have to be wrapped in 
						  const double maximum_bitscore = (std::numeric_limits<double>::max)()) //parenthese sdue to unwanted macro expansion
	{
		auto bitscore_filter = [&](const auto& elem)
		{
			return elem.bitscore >= minimum_bitscore && elem.bitscore <= maximum_bitscore;
		};

		T filtered_results;
		std::ranges::copy(homology_results | std::views::filter(bitscore_filter), std::begin(filtered_results));
		return filtered_results;
	}

	///* Keeps up to @hits_to_keep of best-scoring HMMER results for each query */
	//mipfinder::homology::Results
	//	keepTopHits(const mipfinder::homology::Results& results, std::size_t hits_to_keep);

	///* Removes all entries from HMMER results files where the query is the same
	// * as the target */
	//mipfinder::homology::Results
	//	removeSelfHits(const mipfinder::homology::Results& results);

	///* Returns a list of proteins from @proteins that match the identifiers found
	// * in HMMER @results */
	//mipfinder::ProteinSet
	//	convertToProtein(const mipfinder::homology::Results& results,
	//					 const mipfinder::ProteinSet& proteins);

	////Creates HMMER profiles from all cMIP Multiple Sequence Alignments. 
	////Creates one HMMER profile per one MSA.
	//void createHmmProfiles(const std::filesystem::path& msa_dir,
	//					   const std::filesystem::path& output_dir);

	////Creates one large HMMER profile file from all homologous cMIP HMMER
	////profiles. This simplifies downstream processing.
	//std::filesystem::path
	//	concatenateHmmProfiles(const std::filesystem::path& hmmer_profile_dir,
	//						   const std::filesystem::path& output_file);

	///* Returns all HmmerResults queries where the target length - query length >=
	// * @min_length_difference */
	//mipfinder::homology::Results
	//	filterByLengthDifference(const mipfinder::homology::Results& results,
	//							 const mipfinder::ProteinSet& proteins,
	//							 unsigned int min_length_difference);

	///* Takes HmmerResults and filters out all targets that do not contain the same
	// * protein family identifier as the query */
	//mipfinder::homology::Results
	//	filterByProteinFamilyAndDomain(const mipfinder::homology::Results& results,
	//								   const mipfinder::Proteome& proteome);

	///* Returns all HmmerResults queries that have less or equal than
	// * @maximum_homologues targets */
	//mipfinder::homology::Results
	//	filterByHomologueCount(const mipfinder::homology::Results& results,
	//						   unsigned int maximum_homologues);

}

namespace std
{
	template <>
	struct hash<mipfinder::homology::Result>
	{
		std::size_t operator()(const mipfinder::homology::Result& k) const
		{
			using std::hash;

			return ((hash<std::string>()(k.query)
					^ (hash<std::string>()(k.target) << 1)) >> 1);
		}
	};
}
#endif
