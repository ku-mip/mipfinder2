#ifndef MIPFINDER_MIPFINDER_H
#define MIPFINDER_MIPFINDER_H

#include <filesystem>
#include <string>

#include "configuration.h"
#include "hmmer.h"
#include "hmmer_scoring_data.h"

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

		struct HmmerParameters
		{
			double gap_open_probability;
			double gap_extension_probability;
			std::string scoring_matrix;
		};

		struct FileParamaters
		{
			std::filesystem::path input_proteome;
			std::filesystem::path known_microprotein_list;
			std::filesystem::path interpro_database;
			std::filesystem::path gene_ontology_database;
			std::filesystem::path uniprot_to_intepro_id_conversion_file;
			std::filesystem::path uniprot_to_go_id_conversion_file;
		};

		struct RunParameters
		{
			unsigned int minimum_microprotein_length;
			unsigned int maximum_microprotein_length;
			unsigned int minimum_ancestor_length;
			unsigned int maximum_ancestor_length;
			unsigned int maximum_homologues_per_microprotein;
			unsigned int minimum_length_difference;
			unsigned int maximum_ancestor_count;
			unsigned int maximum_protein_existence_level;
			double microprotein_homologue_bitscore_cutoff;
			double ancestor_bitscore_cutoff;
			std::string output_format;
			std::string organism_identifier;
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
		Configuration configuration;



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
		requires std::ranges::range<T>&& requires (typename std::ranges::range_value_t<T> t)
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
			detail::compareMicroproteinsToMicroproteins(potential_microproteins, hmmer_params, homology_search_output);
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

		//struct FolderParameters
		//{
		//	std::filesystem::path results_folder;
		//	std::filesystem::path msa_folder;
		//	std::filesystem::path hmmprofile_folder;
		//	std::filesystem::path homologue_folder;
		//};

		//HmmerParameters m_hmmer_parameters;
		//FileParamaters m_file_parameters;
		//RunParameters m_run_parameters;

		/* All folders that mipfinder needs to run properly */
		//std::filesystem::path m_results_folder;

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
