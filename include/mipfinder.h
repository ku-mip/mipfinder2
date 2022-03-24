#ifndef MIPFINDER_MIPFINDER_H
#define MIPFINDER_MIPFINDER_H

#include <filesystem>
#include <string>

#include "configuration.h"
#include "hmmer_scoring_data.h"

namespace mipfinder
{
	class Mipfinder
	{
	public:
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

		Mipfinder(const std::filesystem::path& configuration_file);

		//It is responsible for running the actual mipfinder pipeline. See 
		//documentation for details. 
		void run();

		//void writeOutput(std::string filename);

	private:
		enum class Flags
		{
			KNOWN_MIPS_PROVIDED = 0,
			INTERPRO_DATABASE_FILE,
			UNIPROT_TO_INTERPRO_CONVERSION_FILE,
			GO_DATABASE_FILE,
			UNIPROT_TO_GO_CONVERSION_FILE
		};

		Flags run_flags;


		struct ResultsParameters
		{
		};



		struct FolderParameters
		{
			std::filesystem::path results_folder;
			std::filesystem::path msa_folder;
			std::filesystem::path hmmprofile_folder;
			std::filesystem::path homologue_folder;
		};

		HmmerParameters m_hmmer_parameters;
		FileParamaters m_file_parameters;
		RunParameters m_run_parameters;


		/* All folders that mipfinder needs to run properly */
		std::filesystem::path m_results_folder;
		//std::filesystem::path msa_folder_;
		//std::filesystem::path hmmprofile_folder_;
		//std::filesystem::path homologue_folder_;

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
