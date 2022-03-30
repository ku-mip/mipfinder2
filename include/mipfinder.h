#ifndef MIPFINDER_MIPFINDER_H
#define MIPFINDER_MIPFINDER_H

#include <filesystem>
#include <string>

#include "configuration.h"
#include "hmmer.h"
#include "hmmer_scoring_data.h"

namespace detail
{
	template <typename T, typename U>
	requires std::ranges::range<T>&& std::ranges::range<U>
		auto filterAncestorHomologySearchResults(const T& proteome, const U& microprotein_ancestor_homology_results, const mipfinder::Configuration& configuration)
	{
		//Remove homologous results with bitscores outside the acceptable range
		auto homologues = mipfinder::homology::parseResultsFile(microprotein_ancestor_homology_results);

		const double minimum_allowed_homology_bitscore;
		try {
			minimum_allowed_homology_bitscore = std::stod(configuration.value("HMMER", "ancestor_bitscore_cutoff"));
		}
		catch (std::invalid_argument& e) {
			minimum_allowed_homology_bitscore = 10.0;
			LOG(ERROR) << "Ancestor bitscore cutoff specified in configuration file is not a double.";
			LOG(ERROR) << "Using default value of " + std::to_string(minimum_allowed_homology_bitscore) + " instead.";
		}
		catch (std::out_of_range& e) {
			minimum_allowed_homology_bitscore = 10.0;
			LOG(ERROR) << "Ancestor bitscore cutoff specified in configuration file exceeds maximum allowed value.";
			LOG(ERROR) << "Using default value of " + std::to_string(minimum_allowed_homology_bitscore) + " instead.";
		}

		constexpr double maximum_allowed_homology_bitscore = 120;
		auto high_confidence_ancestors = mipfinder::homology::filterByBitscore(homologues, minimum_allowed_homology_bitscore, maximum_allowed_homology_bitscore);

		//Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
		//protein with a similar domain. This would be a problem for very common
		//domains such as kinases, zinc fingers etc.
		constexpr auto& maximum_homologous_ancestors_per_microprotein_to_keep = configuration.value("MIP", "maximum_ancestors");
		auto top_homologues_only = mipfinder::homology::keepTopHomologues(high_confidence_ancestors, maximum_homologous_ancestors_per_microprotein_to_keep);

		//Filter out ancestors that are within x a.a of the cMIP. If a protein is chosen as an ancestor
		//but is only slightly larger than the microprotein, then the ancestor is highly unlikely
		//to contain another domain and cannot function as an ancestor of a microProtein because it
		//does not have another effector domain.
		constexpr auto& minimum_length_difference = configuration.value("MIP", "min_length_difference");

		//if (proteome)
		//filtered_homology_results
	// 
	// 
	// 
	//const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
	//results = mipfinder::homology::filterByLengthDifference(results,
	//														proteome.data(),
	//														min_length_diff);
		return high_confidence_ancestors;
	}

	//Create HMMER profiles for each protein in "homologous_microproteins" that has a homologue.
	template <typename T, typename U>
	void createHmmprofiles(const T& homologous_microproteins, const U& homology_relationship_table, const std::filesystem::path& hmmprofile_output_folder)
	{
		if (!std::filesystem::exists(hmmprofile_output_folder))
		{
			std::filesystem::create_directories(hmmprofile_output_folder);
		}

		LOG(INFO) << "Creating HMMER profiles from homologous microProtein sequences";
		std::size_t profiles_built_counter = 0;
		//Group all unique homologues together into one sequence based on homology
		//for (const auto& [protein_identifier, homologue_identifiers] : homology_relationship_table) {
		//    std::unordered_set<std::ranges::range_value_t<std::remove_const_t<T>>> grouped_homologous_proteins;
		//    if (auto found_protein = std::ranges::find_if(homologous_microproteins, [&protein_identifier](const auto& protein)
		//    {
		//        return protein.identifier() == protein_identifier;
		//    }); found_protein != std::ranges::end(homologous_microproteins)) {
		//        grouped_homologous_proteins.insert(*found_protein);
		//    }
		//    else {
		//        continue;
		//    };

		//    for (const auto& homologue_identifier : homologue_identifiers) {
		//        if (auto found_protein = std::ranges::find_if(homologous_microproteins, [&homologue_identifier](const auto& protein)
		//        {
		//            return protein.identifier() == homologue_identifier.target;
		//        }); found_protein != std::ranges::end(homologous_microproteins)) {
		//            grouped_homologous_proteins.insert(*found_protein);
		//        }
		//    }

		//    //Create files of unaligned sequences as input to an alignment program
		//    const std::filesystem::path unaligned_sequence_dir = hmmprofile_output_folder / "unaligned";
		//    std::filesystem::create_directory(unaligned_sequence_dir);
		//    const std::filesystem::path unaligned_grouped_sequences_fasta = unaligned_sequence_dir / (protein_identifier + ".fasta");
		//    mipfinder::proteinToFasta(grouped_homologous_proteins, unaligned_grouped_sequences_fasta);

		//    //Align the grouped sequences
		//    const std::filesystem::path multiple_sequence_alignment_dir = hmmprofile_output_folder / "msa";
		//    std::filesystem::create_directory(multiple_sequence_alignment_dir);
		//    const std::filesystem::path output_msa_file = multiple_sequence_alignment_dir / (protein_identifier + ".msa");
		//    detail::createMultipleSequenceAlignment(unaligned_grouped_sequences_fasta, output_msa_file);

		//    //Create a hmmprofile from the aligned sequences
		//    const auto hmmprofile_output_file = hmmprofile_output_folder / (protein_identifier + ".hmmprofile");

		//    const auto hmmprofile_name_command = "-n " + protein_identifier; //Names the MSA profile as the protein identifier
		//    mipfinder::homology::buildHmmerProfile(output_msa_file, hmmprofile_output_file, hmmprofile_name_command);
		//    ++profiles_built_counter;
		//}
		LOG(INFO) << "Done creating HMMER profiles";
		LOG(INFO) << "Created " << profiles_built_counter << " profiles";
	}



	template <typename T, typename U>
	void filterOutTooShortAncestors(const T& proteome, const U& homology_search_results, const std::size_t minimum_length_difference)
	{
		U filtered_homology_results{};
		for (const auto& homology_hit : homology_search_results) {
			auto contains_protein = [&homology_hit](const auto& protein)
			{
				return (protein.identifier() == homology_hit.query && protein.identifier() == homology_hit.target);
			};

			if (!std::ranges::find_if(proteome, contains_protein) != std::ranges::end(proteome)) {
				continue;
			}
		}
	}

	//Find which proteins are potential ancestors of the given single-copy microproteins, i.e.
	//which microproteins may have evolved from the larger proteins.
	//
	//@Params
	//single_copy_microproteins - A container 
	template <typename T, typename U, typename V>
	requires std::ranges::range<T>&& std::ranges::range<U>
		void findAncestorsOfSingleCopyMicroproteins(const T& single_copy_microproteins,
			const U& potential_ancestors,
			const V& homology_search_parameters,
			const std::filesystem::path& homology_search_output_file)
	{
		LOG(INFO) << "Finding ancestors of single-copy microProteins";
		if (std::ranges::size(single_copy_microproteins) == 0) {
			LOG(INFO) << "No single-copy microProteins found, aborting homology search";
			return;
		}

		if (std::ranges::size(potential_ancestors) == 0) {
			LOG(INFO) << "No ancestors found, aborting homology search";
			return;
		}

		const std::filesystem::path results_path = homology_search_output_file.parent_path();
		const std::filesystem::path query_file_location = results_path / "single_copy_microproteins.fasta";
		const std::filesystem::path database_location = results_path / "ancestors.fasta";
		mipfinder::protein::proteinToFasta(single_copy_microproteins, query_file_location);
		mipfinder::protein::proteinToFasta(single_copy_microproteins, database_location);

		homology_search_parameters.push_back("--tblout " + homology_search_output_file.string());
		//const auto extra_param = "--mx " + parameters.scoring_matrix;
		//const std::string extra_phmmer_parameters = "--popen " + std::to_string(parameters.gap_open_probability)
		//    + " --pextend " + std::to_string(parameters.gap_extension_probability)
		//    + extra_param;
		mipfinder::homology::phmmer(query_file_location, database_location, homology_search_parameters);
	}

	template <typename T, typename U, typename V>
	requires std::ranges::range<T>&& std::ranges::range<U>
		void findAncestorsOfHomologousMicroproteins(const T& homologous_microproteins,
			const U& potential_ancestors,
			const V& homology_relationship_table,
			const std::filesystem::path& homology_search_output_file)
	{
		const std::filesystem::path parent_folder = homology_search_output_file.parent_path();
		const std::filesystem::path hmmer_output_folder = parent_folder / "hmmprofile";

		detail::createHmmprofiles(homologous_microproteins, homology_relationship_table, hmmer_output_folder);

		LOG(INFO) << "Merging HMMprofiles";
		const std::filesystem::path merged_profile_file = hmmer_output_folder / "mpf_merged.hmmprofile";
		detail::mergeHmmprofileFiles(hmmer_output_folder, merged_profile_file);

		//WIP: Find a better location for this
		const std::filesystem::path ancestor_fasta_file = hmmer_output_folder / "all_ancestors";
		mipfinder::protein::proteinToFasta(potential_ancestors, ancestor_fasta_file);
		LOG(INFO) << "Performing hmmsearch";

		std::filesystem::path database_file{ "hmmsearch_database.txt" };
		std::filesystem::path database_file_location = results_path / database_file;

		std::initializer_list<std::string> options = { "-o /dev/null", "--tblout " + database_file_location.string() + ".txt" };
		mipfinder::homology::hmmsearch(merged_profile_file, ancestor_fasta_file, options);
	}

	template <typename T>
	mipfinder::homology::Results
		filterByLengthDifference(const mipfinder::homology::Results& homology_results,
			const T& proteome,
			const std::size_t min_length_difference)
	{
		mipfinder::homology::Results filtered_results;

		for (const auto& [protein, homologues] : homology_results)
		{
			for (const auto& homologue : homologues) {

			}
		}

		//mipfinder::homology::Results filtered;

		//std::unordered_map<std::string, mipfinder::Protein*> lookup_table;
		//for (const auto& protein : proteins) {
		   // lookup_table[protein->identifier()] = protein;
		//}

		//for (const auto& result : results) {
		   // const auto query = lookup_table.at(result.query);
		   // const auto target = lookup_table.at(result.target);

		   // if (target->length() <= query->length()) {
			  //  continue;
		   // }
		   // assert(target->length() > query->length());

		   // const auto length_difference = target->length() - query->length();
		   // if (length_difference < min_length_difference) {
			  //  continue;
		   // }

		   // filtered.push_back(result);
		//}
		//return filtered;
	}
}

namespace mipfinder
{
	class Mipfinder
	{
	public:
		struct ClassifiedMicroproteins
		{
			mipfinder::protein::ProteinList single_copy;
			mipfinder::protein::ProteinList homologous;
			mipfinder::homology::Results homology_table;
		};

		Mipfinder() = delete;
		Mipfinder(const std::filesystem::path& configuration_file);

		/**
		 *  @brief Start the mipfinder pipeline.
		 */
		void run();

		//void writeOutput(std::string filename);

	private:
		struct HomologyParameters
		{
			double homologue_bitscore_cutoff;
			double ancestor_bitscore_cutoff;
			double gap_open_probability;
			double gap_extend_probability;
			std::string matrix;
		};

		struct MicroproteinParameters
		{
			std::size_t maximum_microprotein_length;
			std::size_t minimum_ancestor_length;
			std::size_t maximum_ancestor_length;
			std::size_t minimum_length_difference;
			std::size_t maximum_ancestor_homologues;
			uint8_t maximum_allowed_protein_existence;
		};

		struct GeneralParameters
		{
			std::string organism_identifier;
		};

		struct FileParameters
		{

		};

		std::filesystem::path configuration_file;
		Configuration configuration;

		HomologyParameters m_hmmer_parameters;
		FileParamaters m_file_parameters;
		RunParameters m_run_parameters;


		//Find all potential microproteins from a proteome based on predetermined criteria
		//
		//@Params
		//proteome - A container of proteins to find microproteins from.
		//run_params - Parameters specifying microprotein characteristics.
		//hmmer_params - Parameters for running homology search.
		//homology_search_output - Path to a folder where the homology search result output
		//                         file will be created.
		//
		//@Return - A struct containing single-copy and homologous microproteins, as well as a homology relationship table
		//for the homologous microproteins.
		template <typename T>
		requires std::ranges::range<T> && requires (typename std::ranges::range_value_t<T> t)
		{
			{ t.sequence().length() } -> std::convertible_to<std::size_t>;
			{ t.identifier().to_string() } -> std::convertible_to<std::string>;
		}
		ClassifiedMicroproteins findMicroproteins(const T& proteins,
												  const std::filesystem::path& homology_search_output)
		{
			LOG(DEBUG) << "Finding microProteins that meet the length criteria as defined in configuration";
			const std::size_t minimum_allowed_microprotein_length = configuration.value("MIP", "minimum_microprotein_length");
			const std::size_t maximum_allowed_microprotein_length = configuration.value("MIP", "maximum_microprotein_length");
			auto potential_microproteins = detail::filterProteinsByLength(
				proteins,
				minimum_allowed_microprotein_length,
				maximum_allowed_microprotein_length);

			if (std::ranges::size(potential_microproteins) == 0) {
				throw std::runtime_error("After filtering the proteome by length, no potential microProteins were found in the proteome. Stopping mipfinder.");
			}

			//Find homologous microproteins
			std::vector<std::string> homology_search_parameters = "--mx " + configuration.value("HMMER", "matrix");
			detail::compareMicroproteinsToMicroproteins(potential_microproteins, homology_search_output);
			//Filter out all microprotein homology results below bitscore_cutoff as these do not denote real
			//homologous relationships
			auto microprotein_homology_results = mipfinder::homology::parseResultsFile(homology_search_output);
			const double lowest_allowed_microprotein_homology_bitscore = configuration.value("HMMER", "homologue_bitscore_cutoff");
			auto strong_homologous_matches = mipfinder::homology::filterByBitscore(microprotein_homology_results, lowest_allowed_microprotein_homology_bitscore);

			//Filter out microproteins with more than @maximum_homologues_allowed homologues. This
			//ensures that we do not pick up large protein families that may contribute to false-positives
			//due to these microproteins containing domains that are very common, e.g. zinc fingers
			const auto& maximum_allowed_homologues_per_microprotein = configuration.value("MIP", "maximum_homologues");
			auto no_large_protein_families = mipfinder::homology::keepTopHomologues(strong_homologous_matches, maximum_allowed_homologues_per_microprotein);

			LOG(DEBUG) << "Exiting filterProteinsByLength()";
			return detail::classifyMicroproteins(potential_microproteins, no_large_protein_families);
		}

		//Creates all necessary results folder to run mipfinder v2.0
		void createFolders();

		//Creates a FASTA file of all homologous cMIP sequences to be used as an 
		//input for clustalo. For each protein it creates a file called 
		//`PROTEIN_ID`_homologues.fasta in `ORGANISM_NAME`/msa/ folder, where
		//`ORGANISM_NAME` is specified in the configuration and `PROTEIN_ID` is the
		//protein UniProt identifier.
		void writeHomologuesToFasta(const mipfinder::protein::ProteinList& homologous_cmips);

		//Runs clustalo on each file in the specified directory.
		void alignHomologousCmips(const std::filesystem::path& unaligned_seq_dir);

		//Performs a hmmsearch on the proteome using a file containing one or more
		//hmmprofiles. Returns a path to the hmmsearch results file which is in a 
		//table format
		std::filesystem::path
			hmmsearchHomologousCmips(const std::filesystem::path& hmmprofiles_file);

		typedef double (*HmmerScoringAlgorithm)(const HmmerScoringData&);

		//void assignHmmerScores(const mipfinder::ProteinSet& proteome,
		//					   HmmerScoringAlgorithm algorithm);
	};
}

#endif
