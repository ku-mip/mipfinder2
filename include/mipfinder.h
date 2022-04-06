#ifndef MIPFINDER_MIPFINDER_H
#define MIPFINDER_MIPFINDER_H

#include <filesystem>
#include <string>

#include "configuration.h"
#include "hmmer.h"
#include "hmmer_scoring_data.h"

namespace detail
{
	//Create HMMER profiles for each protein in "homologous_microproteins" that has a homologue.
	template <typename T, typename U>
	void createHmmprofiles(const T& homologous_microproteins, const U& homology_relationship_table, const std::filesystem::path& hmmprofile_output_folder)
	{
		static_assert(false);
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

	////Find which proteins are potential ancestors of the given single-copy microproteins, i.e.
	////which microproteins may have evolved from the larger proteins.
	////
	////@Params
	////single_copy_microproteins - A container 
	//template <typename T, typename U, typename V>
	//requires std::ranges::range<T>&& std::ranges::range<U>
	//	void findAncestorsOfSingleCopyMicroproteins(const T& single_copy_microproteins,
	//		const U& potential_ancestors,
	//		const V& homology_search_parameters,
	//		const std::filesystem::path& path_to_results_folder)
	//{
	//	LOG(INFO) << "Finding ancestors of single-copy microProteins";
	//	if (std::ranges::size(single_copy_microproteins) == 0) {
	//		LOG(INFO) << "No single-copy microProteins found, stopping homology search";
	//		return;
	//	}

	//	if (std::ranges::size(potential_ancestors) == 0) {
	//		LOG(INFO) << "No ancestors found, stopping homology search";
	//		return;
	//	}

	//	const std::filesystem::path query_file_location = path_to_results_folder / "single_copy_microproteins.fasta";
	//	const std::filesystem::path database_location = path_to_results_folder / "ancestors.fasta";
	//	mipfinder::protein::proteinToFasta(single_copy_microproteins, query_file_location);
	//	mipfinder::protein::proteinToFasta(single_copy_microproteins, database_location);

	//	homology_search_parameters.push_back(" --tblout " + homology_search_output_file.string());
	//	mipfinder::homology::phmmer(query_file_location, database_location, homology_search_parameters);
	//}

	//template <typename T, typename U>
	//	void findAncestorsOfHomologousMicroproteins(const T& homologous_microproteins,
	//		const T& potential_ancestors,
	//		const U& homology_relationship_table,
	//		const std::filesystem::path& homology_search_output_file)
	//{
	//	const std::filesystem::path parent_folder = homology_search_output_file.parent_path();
	//	const std::filesystem::path hmmer_output_folder = parent_folder / "hmmprofile";

	//	detail::createHmmprofiles(homologous_microproteins, homology_relationship_table, hmmer_output_folder);

	//	LOG(INFO) << "Merging HMMprofiles";
	//	const std::filesystem::path merged_profile_file = hmmer_output_folder / "mpf_merged.hmmprofile";
	//	detail::mergeHmmprofileFiles(hmmer_output_folder, merged_profile_file);

	//	//WIP: Find a better location for this
	//	const std::filesystem::path ancestor_fasta_file = hmmer_output_folder / "all_ancestors";
	//	mipfinder::protein::proteinToFasta(potential_ancestors, ancestor_fasta_file);
	//	LOG(INFO) << "Performing hmmsearch";

	//	std::filesystem::path database_file{ "hmmsearch_database.txt" };
	//	std::filesystem::path database_file_location = results_folder / database_file;

	//	std::initializer_list<std::string> options = { "-o /dev/null", "--tblout " + database_file_location.string() + ".txt" };
	//	mipfinder::homology::hmmsearch(merged_profile_file, ancestor_fasta_file, options);
	//}

	//const mipfinder::protein::ProteinList
	//	filterByLengthDifference(const mipfinder::homology::Results& homology_results,
	//		const mipfinder::protein::ProteinList& proteins,
	//		const std::size_t min_length_difference)
	//{
	//	std::unordered_map<std::string, std::vector<std::string>> homology_relation;
	//	for (const auto& result : homology_results)
	//	{
	//		homology_relation[result.query].push_back(result.target);
	//	}

	//	const mipfinder::protein::ProteinList filtered;







	//	std::unordered_map<std::string, mipfinder::protein::Protein*> lookup_table;
	//	for (const auto& protein : proteins) {
	//	    lookup_table[protein->identifier()] = protein;
	//	}

	//	for (const auto& result : homology_relation) {
	//	    const auto query = lookup_table.at(result.query);
	//	    const auto target = lookup_table.at(result.target);

	//	    if (target->length() <= query->length()) {
	//		    continue;
	//	    }

	//	    const auto length_difference = target->length() - query->length();
	//	    if (length_difference < min_length_difference) {
	//		    continue;
	//	    }

	//	    filtered.push_back(result);
	//	}
	//	return filtered;
	//}
}


namespace detail
{
	template <typename T>
	T filterByExistenceLevel(T proteins, const std::size_t max_existence_level)
	{
		auto isAboveAllowedExistenceLevel = [&](const auto& elem) {
			return (elem.existence_level() > max_existence_level);
		};
		std::erase_if(proteins, isAboveAllowedExistenceLevel);
		return proteins;
	}

	template <typename T>
	T filterByLength(T proteins, const std::size_t min_length, const std::size_t max_length)
	{
		auto isNotWithinLengthRange = [&](const auto& elem) {
			return (min_length <= elem.length() <= max_length);
		};
		std::erase_if(proteins, isNotWithinLengthRange);
		return proteins;
	}
}


namespace mipfinder
{
	class Mipfinder
	{
	public:
		struct ClassifiedMicroproteins
		{
			protein::ProteinList single_copy;
			protein::ProteinList homologous;
			homology::Results homology_table;
		};

		struct HomologyParameters
		{
			double homologue_bitscore_cutoff;
			double minimum_ancestor_bitscore;
			double maximum_ancestor_bitscore;
			double gap_open_probability;
			double gap_extend_probability;
			std::string matrix;
		};

		struct MicroproteinParameters
		{
			unsigned int minimum_ancestor_length;
			unsigned int maximum_ancestor_length;
			unsigned int minimum_length_difference;
			unsigned int maximum_ancestor_homologues;
			unsigned int minimum_domains_per_microprotein;
			unsigned int maximum_domains_per_microprotein;
			unsigned int minimum_microprotein_length;
			unsigned int maximum_microprotein_length;
			//unsigned int maximum_homologues_per_microprotein;
			std::optional<unsigned int> maximum_allowed_protein_existence;
		};

		struct GeneralParameters
		{
			std::string organism_identifier;
			std::filesystem::path results_folder;
		};

		struct FileParameters
		{
			std::filesystem::path proteome;
			std::optional<std::filesystem::path> interpro_database;
			std::optional<std::filesystem::path> uniprot_to_intepro;
			std::optional<std::filesystem::path> go_database;
			std::optional<std::filesystem::path> uniprot_to_go;
			std::optional<std::filesystem::path> known_microprotein_identifiers;
		};

		Mipfinder() = delete;
		Mipfinder(const std::filesystem::path& configuration_file);

		/**
		 *  @brief Start the mipfinder pipeline.
		 */
		void run();

		//void writeOutput(std::string filename);
	private:
		std::filesystem::path configuration_file;
		std::filesystem::path results_folder;

		HomologyParameters homology_parameters;
		FileParameters file_parameters;
		GeneralParameters general_parameters;
		MicroproteinParameters microprotein_parameters;

		void setConfiguration(const std::filesystem::path& configuration_file);


		template <typename T>
		T findCandidateMicroproteins(T proteome)
		{
			T candidate_microproteins;
			if (microprotein_parameters.maximum_allowed_protein_existence) {
				LOG(INFO) << "Filtering proteins based on their existence level";
				candidate_microproteins = filterByExistenceLevel(proteome, microprotein_parameters.maximum_allowed_protein_existence.value());
			}
			else {
				LOG(ERROR) << "Maximum protein existence level in configuration is not set (correctly). ";
				LOG(ERROR) << "Skipping filtering based on protein existence level.";
			}

			candidate_microproteins = test::filterByLength(candidate_microproteins,
				microprotein_parameters.minimum_microprotein_length,
				microprotein_parameters.maximum_microprotein_length);

			/* INTERPRO annotation is optional for the analysis */
			if (file_parameters.interpro_database && file_parameters.uniprot_to_intepro) {
				if (!std::filesystem::is_regular_file(file_parameters.uniprot_to_intepro.value())
					|| !std::filesystem::is_regular_file(file_parameters.interpro_database.value())) {
					LOG(ERROR) << "InterPro database or UniProt to InterPro conversion file could not be found";
					LOG(ERROR) << "Skipping InterPro annotation";
				}
				else {
					LOG(INFO) << "Incorporating InterPro data into mipfinder analysis";
					auto interpro_database = mipfinder::Interpro(file_parameters.interpro_database.value());
					auto uniprot_to_interpro_conversion_table = detail::parseProteinDomainList(file_parameters.uniprot_to_intepro.value());

					candidate_microproteins = detail::filterByDomainCount(candidate_microproteins,
						microprotein_parameters.minimum_domains_per_microprotein,
						microprotein_parameters.maximum_domains_per_microprotein,
						interpro_database,
						uniprot_to_interpro_conversion_table);
				}
			}
			return candidate_microproteins;
		}


		template <typename T>
		T findCandidateAncestors(T proteome)
		{
			auto candidate_ancestors = filterByLength(proteome,
				microprotein_parameters.minimum_ancestor_length,
				microprotein_parameters.maximum_ancestor_length);
			return candidate_ancestors;
		}

		/* Classify microProteins into single - copy and homologous microproteins based on the homology search results.
		 * Single-copy microproteins are proteins that do not have any homologues among other microproteins, while
		 * homologous microproteins are proteins that have at least one other homologue among microproteins.
		 *
		 * @Params
		 * microproteins - A container of microproteins to be classified.
		 * homology_search_results - A container of homology search results corresponding to the given @microproteins.
		 *
		 * @Return - Two lists corresponding to single-copy and homologous microproteins. If the @homology_search_results
		 * were not obtained from comparing the @microproteins, the result is undefined. 
		 */
		template <typename T>
		ClassifiedMicroproteins classifyCandidateMicroproteins(const T& candidate_microproteins)
		{
			/* Compare all microproteins to each other to find their evolutionary relationship */
			auto candidate_microproteins_fasta = results_folder / "candidate_microproteins.fasta";
			mipfinder::protein::proteinToFasta(candidate_microproteins, candidate_microproteins_fasta);

			auto homology_search_results = results_folder / "cmips_vs_cmips.txt";
			std::initializer_list<std::string> homology_search_options
			{
				"--mx " + homology_parameters.matrix,
				"--popen " + std::to_string(homology_parameters.gap_open_probability),
				"--pextend " + std::to_string(homology_parameters.gap_extend_probability),
				"-T " + std::to_string(homology_parameters.homologue_bitscore_cutoff),
				"--tblout " + homology_search_results.string()
			};
			mipfinder::homology::phmmer(candidate_microproteins_fasta, candidate_microproteins_fasta, homology_search_options);

			/* Divide the supplied candidate microproteins into single-copy and homologous
			 * microproteins based on the homology search results.
			 */
			std::unordered_map<mipfinder::protein::Identifier, std::size_t> homologue_count_table;
			for (const auto& result : homology_search_results) {
			    ++homologue_count_table[mipfinder::protein::Identifier{result.query}];
			}

			T single_copy_microproteins;
			T homologous_microproteins;
			for (const auto& protein : microproteins) {
				if (!homologue_count_table.contains(protein.identifier())) {
					continue;
				}

				if (homologue_count_table.at(protein.identifier()) == 1) {
					single_copy_microproteins.push_back(protein);
				}
				else {
					homologous_microproteins.push_back(protein);
				}
			}

			return mipfinder::Mipfinder::ClassifiedMicroproteins{ .single_copy = single_copy_microproteins,
													  .homologous = homologous_microproteins,
													  .homology_table = homology_search_results };
		}

		template <typename T>
		mipfinder::homology::Results findAncestorsOfSingleCopyMicroproteins(const T& single_copy_microproteins, const T& ancestors)
		{
			auto single_copy_microproteins_fasta = results_folder / "single_copy_microproteins.fasta";
			mipfinder::protein::proteinToFasta(single_copy_microproteins, single_copy_microproteins_fasta);
			auto ancestors_fasta = results_folder / "ancestors.fasta";
			mipfinder::protein::proteinToFasta(ancestors, ancestors_fasta);

			std::filesystem::path single_vs_ancestors = results_folder / "single_copy_vs_ancestors.txt";
			std::initializer_list<std::string> homology_search_parameters = {
				"--pextend " + std::to_string(homology_parameters.gap_extend_probability),
				"--popen " + std::to_string(homology_parameters.gap_open_probability),
				"--mx " + homology_parameters.matrix,
				"--tblout " + single_vs_ancestors.string()
			};
			mipfinder::homology::phmmer(single_copy_microproteins_fasta, ancestors_fasta, homology_search_parameters);

			/* Results with too low or too high bitscore need to be filtered out. This ensures that 
			 * not too many ancestors with similar domains are picked up by the homology search. 
			 * This can be a problem if the microprotein has a domain that is overrepresented in 
			 * a genome, e.g. specific protein kinase domains.
			 */
			auto filtered_results = mipfinder::homology::parseResultsFile(single_vs_ancestors);
			filtered_results = mipfinder::homology::filterByBitscore(filtered_results,
				homology_parameters.minimum_ancestor_bitscore,
				homology_parameters.maximum_ancestor_bitscore);

			/* If the microprotein has a true ancestor in the proteome, based on their evolution it
			 * is most likely to be one of the top homologues. This step also reduces the amount of
			 * ancestors that may be picked up due to similarity of the microprotein domains to the
			 * larger proteins and simplifies further analysis.
			 */
			mipfinder::homology::keepTopHomologues(filtered_results, microprotein_parameters.maximum_ancestor_homologues);

			return filtered_results;
		}

		template <typename T>
		mipfinder::homology::Results findAncestorsOfHomologousMicroproteins(const T& homologous_microproteins,
			const T& ancestors,
			const mipfinder::homology::Results& microprotein_homology_results)
		{
			const std::filesystem::path hmmer_output_folder = results_folder / "hmmprofile";

			detail::createHmmprofiles(homologous_microproteins, microprotein_homology_results, hmmer_output_folder);

			LOG(INFO) << "Merging HMMprofiles";
			const std::filesystem::path merged_profile_file = hmmer_output_folder / "merged_profiles.hmmprofile";
			detail::mergeHmmprofileFiles(hmmer_output_folder, merged_profile_file);

			const std::filesystem::path ancestor_fasta_file = hmmer_output_folder / "ancestors.fasta";
			mipfinder::protein::proteinToFasta(ancestors, ancestor_fasta_file);
			LOG(INFO) << "Performing hmmsearch";

			std::filesystem::path homologous_vs_ancestors_results = results_folder / "hmmsearch_results.txt";

			std::initializer_list<std::string> options{
				"-o /dev/null",
				"--tblout " + homologous_vs_ancestors_results.string()
			};
			mipfinder::homology::hmmsearch(merged_profile_file, ancestor_fasta_file, options);

			auto filtered_results = mipfinder::homology::parseResultsFile(homologous_vs_ancestors_results);
			filtered_results = mipfinder::homology::filterByBitscore(homologous_vs_ancestors_results,
				homology_parameters.minimum_ancestor_bitscore,
				homology_parameters.maximum_ancestor_bitscore);

			mipfinder::homology::keepTopHomologues(filtered_results, microprotein_parameters.maximum_ancestor_homologues);

			return filtered_results;
		}


		//Creates all necessary results folder to run mipfinder v2.0
		void createFolders();

		//mipfinder::protein::ProteinList
		//filterAncestorHomologySearchResults(const mipfinder::protein::ProteinList& proteome,
		//	const std::filesystem::path& microprotein_ancestor_homology_results);




		//Creates a FASTA file of all homologous cMIP sequences to be used as an 
		//input for clustalo. For each protein it creates a file called 
		//`PROTEIN_ID`_homologues.fasta in `ORGANISM_NAME`/msa/ folder, where
		//`ORGANISM_NAME` is specified in the configuration and `PROTEIN_ID` is the
		//protein UniProt identifier.
		void writeHomologuesToFasta(const protein::ProteinList& homologous_cmips);

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
