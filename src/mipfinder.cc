#include <cassert>
#include <concepts>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <set>
#include <unordered_set>

#include "aliases.h"
#include "ancestor.h"
#include "configuration.h"
#include "easylogging++.h"
#include "fasta.h"
#include "file.h"
#include "go.h"
#include "hmmer.h"
#include "helpers.h"
#include "interpro.h"
#include "mipfinder.h"
#include "protein.h"
#include "proteome.h"
#include "printer.h"

#include <ranges>

using ProteinList = std::vector<Protein>;

INITIALIZE_EASYLOGGINGPP

/* Helper functions */
namespace
{

	/* Ensures that all files required by mipfinder are found. If some
	 * dependencies are missing, throws a std::runtime_error */
	void checkFileDependencies(mipfinder::Configuration& config)
	{
		std::filesystem::path proteome_fasta = config["TARGET"]["organism_proteins_fasta"];
		std::filesystem::path go_id_descriptions = config["GO"]["go_database"];
		std::filesystem::path uniprot_id_to_go = config["GO"]["uniprot_to_go"];
		std::filesystem::path interpro_database = config["INTERPRO"]["interpro_database"];
		std::filesystem::path uniprot_id_to_interpro = config["INTERPRO"]["uniprot_to_interpro"];

		mipfinder::file::exists(proteome_fasta);
		mipfinder::file::exists(go_id_descriptions);
		mipfinder::file::exists(uniprot_id_to_go);
		mipfinder::file::exists(interpro_database);
		mipfinder::file::exists(uniprot_id_to_interpro);
	}

	void configureLogger()
	{
		el::Configurations default_logging;
		default_logging.setToDefault();
		default_logging.set(el::Level::Debug,
							el::ConfigurationType::Enabled,
							"true");
		default_logging.set(el::Level::Debug,
							el::ConfigurationType::Filename,
							"debug.log");
		default_logging.set(el::Level::Info,
							el::ConfigurationType::Filename,
							"mipfinder.log");
		default_logging.set(el::Level::Info,
							el::ConfigurationType::Format,
							"%datetime %level %msg");

		el::Loggers::reconfigureLogger("default", default_logging);
	}

	double scoring_algorithm(const mipfinder::HmmerScoringData& data)
	{
		std::size_t cmip_length = data.query->length();
		std::size_t ancestor_length = data.target->length();
		double cmip_coverage = static_cast<double>(cmip_length / ancestor_length);

		double growth_rate = 0.1;
		double midpoint = 50.0;
		double weight = 1 / (1 + std::exp(growth_rate * (cmip_coverage * 100 - midpoint)));

		double assigned_score = (1.0 / std::pow(data.bitscore, cmip_coverage))
			* data.bitscore
			* weight;
		return assigned_score;
	}

	/* For each cMIP in @results (HMMER queries), add the identified ancestor
	 * (HMMER targets) to the cMIP. Does not add self as an ancestor */
	void associateAncestorsWithCmips(const mipfinder::Proteome& proteome,
									 const mipfinder::HmmerResults& results)
	{
		for (const auto& result : results) {
			if (result.query == result.target) {
				continue;
			}

			const auto cmip = proteome.find(result.query);
			const auto ancestor = proteome.find(result.target);

			/* If ancestors do not have any domains or family annotations, ignore them */
			const auto all_ancestor_domains = ancestor->interproEntries();
			unsigned int domain_family_count{0};
			for (const auto& domain : all_ancestor_domains)
				if (domain.type == mipfinder::Interpro::Type::DOMAIN_TYPE ||
					domain.type == mipfinder::Interpro::Type::FAMILY) {
					++domain_family_count;
				}

			if (domain_family_count == 0) {
				continue;
			}

			assert(cmip != nullptr);
			assert(ancestor != nullptr);
			cmip->addAncestor(mipfinder::Ancestor{ancestor, result.bitscore});
		}
	}

	/* For each cMIP in @results (HMMER queries), add the identified homologues
	 * (HMMER targets, proteins) to the cMIP. Does not add self as a homologue */
	void associateHomologuesWithCmips(const mipfinder::ProteinSet& cmips,
									  const mipfinder::HmmerResults& results)
	{
		/* Create a lookup table */
		std::unordered_map<std::string, mipfinder::Protein*> lookup_table;
		for (const auto& cmip : cmips) {
			lookup_table[cmip->identifier()] = cmip;
		}

		for (const auto& result : results) {
			if (result.query == result.target) {
				continue;
			}

			if (lookup_table.count(result.query) == 0) {
				continue;
			}

			if (lookup_table.count(result.target) == 0) {
				continue;
			}

			const auto cmip = lookup_table[result.query];
			const auto homologue = lookup_table[result.target];
			assert(cmip != nullptr);
			assert(homologue != nullptr);
			cmip->addHomologue(mipfinder::Homologue{homologue, result.bitscore});
		}
	}

	mipfinder::HmmerResults
		filterSingleDomainAncestors(const mipfinder::Proteome& proteome,
									const mipfinder::HmmerResults& results)
	{
		mipfinder::HmmerResults filtered;

		for (const auto& result : results) {
			const auto ancestor_id = result.target;
			const auto ancestor = proteome.find(ancestor_id);
			assert(ancestor != nullptr);

			unsigned domain_count{0};
			for (const auto& interpro_entry : ancestor->interproEntries()) {
				if (interpro_entry.type == mipfinder::Interpro::Type::DOMAIN_TYPE ||
					interpro_entry.type == mipfinder::Interpro::Type::REPEAT) {
					++domain_count;
				}
			}

			/* If it is 0, it is very likely that the annotation of that gene has
			* failed since virtually all proteins are made up of small domains.
			* If it is exactly 1, we assume that the annotation is correct and it
			* really only does have one domain */
			if (domain_count != 1) {
				filtered.push_back(result);
			}
		}
		return filtered;
	}

	/* Divides a proteome into potential cMIPs and potential ancestor protein
	 * based on lengths specified in the configuration file. */
	std::pair<mipfinder::ProteinSet, mipfinder::ProteinSet>
		divideProteome(const mipfinder::Proteome& proteome,
					   const mipfinder::Configuration& config)
	{
		std::size_t max_mip_length = std::stoi(config["MIP"]["max_mip_length"]);
		std::size_t min_ancestor_length = std::stoi(config["MIP"]["min_ancestor_length"]);

		const auto cmips = mipfinder::filterByLength(proteome.data(),
													 1,
													 max_mip_length);
		const auto ancestors = mipfinder::filterByLength(proteome.data(),
														 min_ancestor_length);

		/* Mark respective parts of the proteome as a cMIP or an ancestor */
		for (const auto& cmip : cmips) {
			assert(cmip->type() == mipfinder::Protein::Type::UNKNOWN);
			cmip->setType(mipfinder::Protein::Type::CMIP);
		}
		for (const auto& ancestor : ancestors) {
			assert(ancestor->type() == mipfinder::Protein::Type::UNKNOWN);
			ancestor->setType(mipfinder::Protein::Type::ANCESTOR);
		}

		return std::make_pair(cmips, ancestors);
	}

	/* Applies a set of filters to results file that represent cMIPS being
	 * searched against potential ancestors. Returns a list of filtered HMMER
	 * results */
	mipfinder::HmmerResults filterAncestors(mipfinder::HmmerResults results,
											const mipfinder::Proteome& proteome,
											const mipfinder::Configuration& config)
	{
		/* Apply filters to the unique cMIP vs ancestor results to eliminate
		* low-confidence results */
		const double ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]);
		results = mipfinder::hmmer::filterByBitscore(results,
													 ancestor_bitscore_cutoff);

		results = mipfinder::hmmer::filterByBitscore(results,
													 120,
													 std::less_equal<>());

		/* Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
		* protein with a similar domain. This would be a problem for very common
		* domains such as kinases, zinc fingers etc. */
		const std::size_t hits_to_keep = std::stoi(config["MIP"]["max_ancestor_count"]);
		results = mipfinder::hmmer::keepTopHits(results, hits_to_keep);

		//Filter out ancestors that are within 40 a.a of the cMIP
		const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
		results = mipfinder::hmmer::filterByLengthDifference(results,
															 proteome.data(),
															 min_length_diff);

		results = filterSingleDomainAncestors(proteome, results);

		return results;
	}

	void setCmipTypes(const mipfinder::Proteome& proteome)
	{
		for (const auto& protein : proteome.data()) {
			if (protein->type() != mipfinder::Protein::Type::CMIP) {
				continue;
			}
			if (protein->homologues().empty()) {
				protein->setType(mipfinder::Protein::Type::SINGLE_COPY);
			}
			else {
				protein->setType(mipfinder::Protein::Type::HOMOLOGOUS);
			}
		}
	}

}


namespace detail
{
	mipfinder::ProteinList loadProteome(const std::filesystem::path& fasta_file)
	{
		auto file = mipfinder::file::open(fasta_file);
		const mipfinder::FastaRecords proteome_fasta_records = mipfinder::fasta::extractRecords(file);

		mipfinder::ProteinList proteome;

		for (const auto& [header, sequence] : proteome_fasta_records) {

			const auto& [protein_id, sequence_version, description, existence_level] = mipfinder::fasta::extractUniprotHeader(header);
			const std::string identifier = protein_id + "." + sequence_version;

			proteome.emplace_back(mipfinder::Protein{identifier, sequence, description, std::stoi(existence_level)});
		}
		return proteome;
	}

	template <typename T>
	concept ProteinRange = std::ranges::range<T> && std::same_as<typename std::ranges::range_value_t<T>, mipfinder::Protein>;

	template <typename T>
	concept HasIdentifier = requires (T t)
	{
		{ t.identifier() } -> std::convertible_to<std::string>;
	};

	template <typename T>
	concept HasSequence = requires (T t)
	{
		{ t.sequence() } -> std::convertible_to<std::string>;
	};


	//Creates a FASTA file from existing Proteins
	template <typename T>
	requires std::ranges::range<T>&& HasIdentifier<std::ranges::range_value_t<T>>&& HasSequence<std::ranges::range_value_t<T>>
		void proteinToFasta(T& proteins, const std::filesystem::path& output)
	{
		std::ofstream f;
		f.open(output);
		for (const auto& protein : proteins) {
			f << ">" << protein.identifier() << "\n" << protein.sequence() << "\n";
		}
	}

	//mipfinder::ScoringMatrix detectScoringMatrix(std::string matrix)
	//{
	//	for (auto& c : matrix) {
	//		c = std::toupper(c);
	//	}
	//	if (matrix == "BLOSUM50") {
	//		return mipfinder::ScoringMatrix::BLOSUM50;
	//	}
	//	else if (matrix == "BLOSUM62") {
	//		return mipfinder::ScoringMatrix::BLOSUM62;
	//	}
	//	throw std::logic_error("Invalid scoring matrix option detected, aborting.");
	//}

	//Classify all candiate microProteins as either single-copy (no homologous proteins
	//in the proteome) or homologous (has homologues in the proteome).
	template <typename T, typename U>
	requires std::ranges::range<T>&& HasIdentifier<std::ranges::range_value_t<T>>&& HasSequence<std::ranges::range_value_t<T>>&& std::convertible_to<U, std::filesystem::path>
		void classifyMicroproteins(T& candidate_microproteins, mipfinder::HmmerParameters parameters, const U& phmmer_results_output)
	{
		//FOR NOW: Use files to call phmmer, in the future pipe the data in from stdin (to phmmer)
		LOG(INFO) << "Grouping cMIPs based on homology";
		//Create the input file for phmmer in FASTA format
		const std::string microproteins_fasta_file{"all_microproteins.txt"};
		detail::proteinToFasta(candidate_microproteins, microproteins_fasta_file);
		const auto extra_param = " --mx " + parameters.scoring_matrix;
		mipfinder::hmmer::phmmer(microproteins_fasta_file, microproteins_fasta_file, phmmer_results_output, extra_param);
		LOG(INFO) << "Finished grouping cMIPS";

		////Pipe in every single microProtein into phmmer one by one to compare against the database
		//for (const auto& protein : candidate_microproteins)
		//{
		//#if defined _WIN64 
		//    //Not implemented yet
		//    throw std::logic_error("No windows implementation...yet");


		//#elif defined unix || __unix__ || __unix || __linux__


		//    //int pipefd[2];
		//    //pid_t cpid;
		//    //char buf;
		//    //if (argc != 2) {
		//    //    fprintf(stderr, "Usage: %s <string>\n", argv[0]);
		//    //    exit(EXIT_FAILURE);
		//    throw std::logic_error("Implement UNIX piping");

		//#endif
		//}
	}
}



mipfinder::Mipfinder::Mipfinder(const std::filesystem::path& configuration_file)
	: config_(Configuration(configuration_file)),
	proteome_(Proteome(config_["TARGET"]["organism_proteins_fasta"]))
{
	configureLogger();

	LOG(DEBUG) << "Checking file dependenies";
	checkFileDependencies(config_);

	LOG(DEBUG) << "Setting up GO database";
	mipfinder::Go go_database{config_["GO"]["go_database"]};
	LOG(DEBUG) << "Adding GO ids to proteins";
	addGoIdentifiers(proteome_,
					 go_database,
					 config_["GO"]["uniprot_to_go"]);

	LOG(DEBUG) << "Setting up InterPro database";
	mipfinder::Interpro interpro_database{config_["INTERPRO"]["interpro_database"]};
	LOG(DEBUG) << "Adding InterPro ids to proteins";
	addInterproIdentifiers(proteome_,
						   interpro_database,
						   config_["INTERPRO"]["uniprot_to_interpro"]);

	createFolders();

	//Extract run configuration
	Configuration config{configuration_file};

	mipfinder::HmmerParameters m_hmmer_parameters{
		.homologue_bitscore_cutoff = std::stod(config["HMMER"]["matrix"]),
		.ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]),
		.gap_open_probability = std::stod(config["HMMER"]["gap_open_probability"]),
		.gap_extension_probability = std::stod(config["HMMER"]["gap_extend_probability"]),
		.scoring_matrix = config["HMMER"]["matrix"]
	};





		//auto hmmer_parameters = detail::extractHmmerParameters(configuration_file);
		//auto mipfinder_parameters = detail::extractMipfinderParameters(configuration_file);
		//auto result_files_parameters = detail::extractResultFilesParameters(configuration_file);


		/* Provide a copy of the run parameters with the results */
	//std::filesystem::path config_filename = std::filesystem::path{configuration_file}.filename();
	//std::filesystem::path config_file_copy = results_folder_ / config_filename;
	//std::filesystem::copy(configuration_file,
	//					  config_file_copy,
	//					  std::filesystem::copy_options::overwrite_existing);


}



namespace mipfinder
{

	//////////////////////////////
	//                          //
	//                          //
	//      MIPFINDER V2.0      //
	//                          //
	//                          //
	//////////////////////////////

	void Mipfinder::run()
	{
		LOG(INFO) << "Starting mipfinder v2.0";
		LOG(INFO) << "Detected " << proteome_.size() << " proteins";

		auto proteome = detail::loadProteome(config_["TARGET"]["organism_proteins_fasta"]);

		//Filter out proteins whose existence level hints suggests that they are not translated transcripts
		const unsigned int kMaximumAllowedExistenceLevel = std::stoi(config_["TARGET"]["protein_existence"]);
		auto protein_existence_filter = [&](const auto& protein)
		{
			return protein.existenceLevel() <= kMaximumAllowedExistenceLevel;
		};
		auto real_proteins = proteome | std::views::filter(protein_existence_filter);

		//Divide proteome into two sets, microProteins and ancestors, based on their length
		const std::size_t kMaximumMicroproteinLength = std::stoi(config_["MIP"]["max_mip_length"]);
		auto microprotein_filter = [&](const auto& protein) { return protein.length() <= kMaximumMicroproteinLength; };
		auto candidate_microproteins = real_proteins | std::views::filter(microprotein_filter);

		const std::size_t kMinimumAncestorLength = std::stoi(config_["MIP"]["min_ancestor_length"]);
		auto ancestor_filter = [&](const auto& protein) { return protein.length() <= kMinimumAncestorLength; };
		auto candidate_ancestors = real_proteins | std::views::filter(ancestor_filter);


		const std::filesystem::path classified_microproteins = results_folder_ / "all_cmips_vs_cmips.txt";
		detail::classifyMicroproteins(candidate_microproteins, m_hmmer_parameters, classified_microproteins);

		/////WIP UNDERNEATH


		/* Filter out all results below @bitscore_cutoff as these do not denote real
		 * homologous relationships */
		const double kLowestHomologyBitscoreCutoff = std::stod(config_["HMMER"]["homologue_bitscore_cutoff"]);
		auto bitscore_filter = [&](const auto& protein)
		{
			return protein.bitscore() >= kLowestHomologyBitscoreCutoff;
		};








		//----------------------------------------------------


	  /* The following are the steps in mipfinder v2.0 algorithm */
		auto [cmips, ancestors] = divideProteome(proteome_, config_);
		LOG(INFO) << "Found " << cmips.size() << " cMIPS and " << ancestors.size() << " ancestors";

		int max_protein_existence = std::stoi(config_["TARGET"]["protein_existence"]);
		cmips = filterByExistence(cmips, max_protein_existence);
		ancestors = filterByExistence(ancestors, max_protein_existence);

		/* Find single-copy and homologous cMIPS */
		const auto cmip_hmmer_file = phmmerAgainstSelf(cmips);
		auto cmip_hmmer_results = hmmer::parseResultsFile(cmip_hmmer_file);

		/* Filter out all results below @bitscore_cutoff as these do not denote real
		 * homologous relationships */
		double homologue_bitscore_cutoff = std::stod(config_["HMMER"]["homologue_bitscore_cutoff"]);
		auto high_confidence_cmips = hmmer::filterByBitscore(cmip_hmmer_results,
															 homologue_bitscore_cutoff);

		/* Associate homologues from HMMER results with Protein objects */
		associateHomologuesWithCmips(cmips, high_confidence_cmips);

		mipfinder::ProteinSet unique_cmips;
		mipfinder::ProteinSet homologous_cmips;

		/* Find all unique and homologous cMIPS */
		for (const auto& protein : cmips) {
			if (protein->homologues().empty()) {
				unique_cmips.insert(protein);
			}
			else {
				homologous_cmips.insert(protein);
			}
		}

		/* Deal with single copy cMIPS */
		/* Take all single-copy cMIPS and phmmer them against large proteins to find
		 * their ancestors */
		const auto unique_vs_ancestor_file =
			phmmerCmipsVsAncestors(unique_cmips, "unique_vs_proteome.txt");

		auto unique_vs_ancestor_results = hmmer::parseResultsFile(unique_vs_ancestor_file);

		unique_vs_ancestor_results = filterAncestors(unique_vs_ancestor_results,
													 proteome_,
													 config_);

		associateAncestorsWithCmips(proteome_,
									unique_vs_ancestor_results);

		/* Deal with homologous cMIPS */
		/* Build HMM profiles of all homologous cMIPS and hmmsearch these profiles
		 * against the ancestors */
		writeHomologuesToFasta(homologous_cmips);
		alignHomologousCmips(homologue_folder_);

		if (config_["DEBUG"]["create_hmm_profiles"] == "true") {
			hmmer::createHmmProfiles(msa_folder_, hmmprofile_folder_);
		}

		std::filesystem::path filename{"all_hmmprofiles.txt"};
		std::filesystem::path unified_hmmprofile_file = results_folder_ / filename;
		const auto all_hmmer_profiles = hmmer::concatenateHmmProfiles(hmmprofile_folder_,
																	  unified_hmmprofile_file);

		/* Merges all individual hmm profile files into one large one */
		const auto homologous_vs_ancestor_search =
			hmmsearchHomologousCmips(all_hmmer_profiles);

		auto homologous_vs_ancestor_results = hmmer::parseResultsFile(homologous_vs_ancestor_search);

		homologous_vs_ancestor_results = filterAncestors(homologous_vs_ancestor_results,
														 proteome_,
														 config_);

		associateAncestorsWithCmips(proteome_,
									homologous_vs_ancestor_results);

		setCmipTypes(proteome_);

		/* Assign KNOWN_MIP type to all known mips */
		Proteome known_mips(config_["MIP"]["known_mips_fasta"]);
		for (const auto& mip : known_mips.data()) {
			auto mip_in_proteome = proteome_.find(mip->identifier());
			if (mip_in_proteome) {
				mip_in_proteome->setType(mipfinder::Protein::Type::KNOWN_MIP);
			}
		}

		/* Filter out proteins with more than @maximum_homologues_allowed homologues */
		mipfinder::ProteinSet proteins_to_score;
		for (const auto& protein : proteome_.data()) {
			const double max_homologues_allowed = std::stod(config_["MIP"]["maximum_homologues"]);
			if (protein->homologues().size() > max_homologues_allowed) {
				continue;
			}
			proteins_to_score.insert(protein);
		}

		assignHmmerScores(proteins_to_score, &scoring_algorithm);
		assignHmmerScores(proteins_to_score, &scoring_algorithm);

		/* Filter out proteins with a 0 score */
		mipfinder::ProteinSet non_zero;
		for (const auto& protein : proteome_.data()) {
			if (protein->score() != 0) {
				non_zero.insert(protein);
			}
		}

		/* Filter out proteins that are not CMIPS */
		mipfinder::ProteinSet only_cmip_type;
		for (const auto& protein : non_zero) {
			if ((protein->type() == mipfinder::Protein::Type::SINGLE_COPY) ||
				(protein->type() == mipfinder::Protein::Type::HOMOLOGOUS) ||
				(protein->type() == mipfinder::Protein::Type::KNOWN_MIP)) {
				only_cmip_type.insert(protein);
			}
		}

		/* Sort proteome based on score */
		auto comp = [](const mipfinder::Protein* lhs, const mipfinder::Protein* rhs)
		{
			return lhs->score() > rhs->score();
		};
		std::vector<Protein*> sorted_proteome{only_cmip_type.begin(), only_cmip_type.end()};
		std::sort(sorted_proteome.begin(), sorted_proteome.end(), comp);

		std::filesystem::path results_fasta_file{"cmip_results.fasta"};
		std::filesystem::path results_fasta_location = results_folder_ / results_fasta_file;
		mipfinder::proteinToFasta(sorted_proteome, results_fasta_location);

		LOG(INFO) << "Writing final report";
		const std::string report_format = config_["REPORT"]["format"];
		std::filesystem::path final_results = results_folder_ / "final_results.txt";
		mipfinder::printer::createReport(report_format,
										 '\t',
										 sorted_proteome,
										 final_results);

		LOG(INFO) << "mpf v2.0 has finished.";
	}

	void Mipfinder::createFolders()
	{
		//Creates the main results folder for the run
		const std::string organism_id = config_["TARGET"]["organism_identifier"];
		const std::string folder_name = "results_" + organism_id;

		const std::filesystem::path results_folder{folder_name};
		std::filesystem::create_directory(results_folder);
		results_folder_ = results_folder;

		//Create the subfolders for individual software package results
		//createResultsFolder() has to be run before this is called, as it sets `results_folder_`
		msa_folder_ = results_folder_ / std::filesystem::path{"msa"};
		hmmprofile_folder_ = results_folder_ / std::filesystem::path{"hmmprofile"};
		homologue_folder_ = results_folder_ / std::filesystem::path{"homologues"};
		std::filesystem::create_directory(msa_folder_);
		std::filesystem::create_directory(hmmprofile_folder_);
		std::filesystem::create_directory(homologue_folder_);
	}

	std::filesystem::path Mipfinder::phmmerAgainstSelf(const mipfinder::ProteinSet& cmips)
	{
		std::filesystem::path results_file{"all_cmips_against_cmips.txt"};
		const auto hmmer_cmips_vs_cmips = results_folder_ / results_file;

		if (config_["DEBUG"]["phmmer_cmips"] == "false") {
			return std::filesystem::path{hmmer_cmips_vs_cmips};
		}

		LOG(INFO) << "Grouping cMIPs based on homology";

		const auto scoring_matrix = config_["HMMER"]["matrix"];
		const auto extra_param = " --mx " + scoring_matrix;

		mipfinder::hmmer::phmmer(cmips, cmips, hmmer_cmips_vs_cmips.string(), extra_param);
		LOG(INFO) << "Finished grouping cMIPS";
		return hmmer_cmips_vs_cmips;
	}

	std::filesystem::path Mipfinder::phmmerCmipsVsAncestors(const mipfinder::ProteinSet& cmips,
															const std::filesystem::path& results_filename)
	{
		unsigned int min_ancestor_length =
			std::stoi(config_["MIP"]["min_ancestor_length"]);

		const auto ancestors =
			mipfinder::filterByLength(proteome_.data(), min_ancestor_length);

		const auto scoring_matrix = config_["HMMER"]["matrix"];
		const auto extra_param = " --mx " + scoring_matrix;

		const auto gap_open_probability = config_["HMMER"]["gap_open_probability"];
		const auto gap_extend_probability = config_["HMMER"]["gap_extend_probability"];
		const std::string extra_phmmer_parameters = "--popen "
			+ gap_open_probability
			+ " --pextend "
			+ gap_extend_probability
			+ extra_param;

		std::filesystem::path results_file_location =
			results_folder_ / results_filename;

		if (config_["DEBUG"]["find_unique_cmip_ancestors"] == "true"
			&& config_["DEBUG"]["find_homologous_cmip_ancestors"] == "true") {
			mipfinder::hmmer::phmmer(cmips, ancestors, results_file_location.string());
		}
		return results_file_location;
	}

	void Mipfinder::writeHomologuesToFasta(const mipfinder::ProteinSet& homologous_cmips)
	{
		if (config_["DEBUG"]["write_homologous_groups_fasta"] == "false") { return; }

		LOG(INFO) << "Creating FASTA files of homologous cMIPs";
		int grouped_fasta_files{0};
		for (const auto& protein : homologous_cmips) {
			mipfinder::ProteinSet grouped_homologues;
			grouped_homologues.insert(protein);

			for (const auto& homologue : protein->homologues()) {
				grouped_homologues.insert(homologue.protein);
			}

			const auto id = protein->identifier();

			const std::filesystem::path unaligned_filename{id + "_homologues.fasta"};
			const std::filesystem::path unaligned_file_location{homologue_folder_
				/ unaligned_filename};

			proteinToFasta(grouped_homologues, unaligned_file_location.string());
			++grouped_fasta_files;
		}
		LOG(DEBUG) << "Created " << grouped_fasta_files << " homologous cMIP sequence "
			<< "FASTA files";
	}

	void Mipfinder::alignHomologousCmips(const std::filesystem::path& unaligned_seq_folder)
	{
		if (config_["DEBUG"]["align_homologous_cmips"] == "false") { return; }

		LOG(INFO) << "Aligning groups of homologous cMIPs";

		typedef std::filesystem::directory_iterator DirectoryIter;

		int msas_created{0};
		for (const auto& file : DirectoryIter{unaligned_seq_folder}) {

			//Unaligned FASTA files are named as `PROTEIN_ID`_homologues.fasta. We
			//need to extract the `PROTEIN_ID` part for the aligned filename.
			const auto unaligned_filename = file.path().filename().string();
			const auto tokens = mipfinder::tokenise(unaligned_filename, '_');
			const auto protein_id = tokens[0];

			const std::filesystem::path unaligned_file_path = file.path();

			const std::filesystem::path msa_filename{protein_id + "_aligned.fasta"};
			const std::filesystem::path msa_file_full_path =
				msa_folder_ / msa_filename;

			const std::string clustalo_command = "clustalo -i "
				+ unaligned_file_path.string()
				+ " -o "
				+ msa_file_full_path.string();

			int sys_call_result = std::system(clustalo_command.c_str());
			if (sys_call_result != 0) {
				continue;
			}
			++msas_created;
		}
		LOG(INFO) << "Created " << msas_created << " MSA from grouped cMIP homologues";
	}

	std::filesystem::path
		Mipfinder::hmmsearchHomologousCmips(const std::filesystem::path& hmmprofiles_file)
	{
		unsigned int min_ancestor_length =
			std::stoi(config_["MIP"]["min_ancestor_length"]);

		const auto ancestors =
			mipfinder::filterByLength(proteome_.data(), min_ancestor_length);

		const std::filesystem::path results_filename{"homologous_vs_proteome.txt"};
		const std::filesystem::path results_file_location =
			results_folder_ / results_filename;

		if (config_["DEBUG"]["find_homologous_cmip_ancestors"] == "true") {
			mipfinder::hmmer::hmmsearch(hmmprofiles_file, ancestors, results_file_location);
		}
		return results_file_location;
	}

	void Mipfinder::assignHmmerScores(const mipfinder::ProteinSet& proteins,
									  HmmerScoringAlgorithm algorithm)
	{
		int score_counter{0};
		for (const auto& protein : proteins) {
			/* Only cMIPs have ancestors */
			if (protein->ancestors().empty()) {
				continue;
			}

			for (const auto& ancestor : protein->ancestors()) {
				const double bitscore = ancestor.bitscore;

				double score = algorithm(HmmerScoringData{protein, ancestor.protein, bitscore});
				protein->changeScore(score);
				++score_counter;
			}
		}
	}

}
