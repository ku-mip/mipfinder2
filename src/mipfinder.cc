#include <cassert>
#include <concepts>
#include <cmath>
#include <filesystem>
#include <functional>
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
	void checkFileDependencies(mipfinder::Mipfinder::FileParamaters file_parameters)
	{
		mipfinder::file::exists(file_parameters.input_proteome);
		mipfinder::file::exists(file_parameters.known_microprotein_list);
		mipfinder::file::exists(file_parameters.gene_ontology_database);
		mipfinder::file::exists(file_parameters.interpro_database);
		mipfinder::file::exists(file_parameters.uniprot_to_intepro_id_conversion_file);
		mipfinder::file::exists(file_parameters.uniprot_to_go_id_conversion_file);
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
	//void associateAncestorsWithCmips(const mipfinder::Proteome& proteome,
	//								 const mipfinder::homology::Results& results)
	//{
	//	for (const auto& result : results) {
	//		if (result.query == result.target) {
	//			continue;
	//		}

	//		const auto cmip = proteome.find(result.query);
	//		const auto ancestor = proteome.find(result.target);

	//		/* If ancestors do not have any domains or family annotations, ignore them */
	//		const auto all_ancestor_domains = ancestor->interproEntries();
	//		unsigned int domain_family_count{0};
	//		for (const auto& domain : all_ancestor_domains)
	//			if (domain.type == mipfinder::Interpro::Type::DOMAIN_TYPE ||
	//				domain.type == mipfinder::Interpro::Type::FAMILY) {
	//				++domain_family_count;
	//			}

	//		if (domain_family_count == 0) {
	//			continue;
	//		}

	//		assert(cmip != nullptr);
	//		assert(ancestor != nullptr);
	//		cmip->addAncestor(mipfinder::Ancestor{ancestor, result.bitscore});
	//	}
	//}

	///* For each cMIP in @results (HMMER queries), add the identified homologues
	// * (HMMER targets, proteins) to the cMIP. Does not add self as a homologue */
	//void associateHomologuesWithCmips(const mipfinder::ProteinSet& cmips,
	//								  const mipfinder::homology::Results& results)
	//{
	//	/* Create a lookup table */
	//	std::unordered_map<std::string, mipfinder::Protein*> protein_lookup_table;
	//	for (const auto& cmip : cmips) {
	//		protein_lookup_table[cmip->identifier()] = cmip;
	//	}

	//	for (const auto& result : results) {
	//		if (result.query == result.target) {
	//			continue;
	//		}

	//		if (protein_lookup_table.count(result.query) == 0) {
	//			continue;
	//		}

	//		if (protein_lookup_table.count(result.target) == 0) {
	//			continue;
	//		}

	//		const auto cmip = protein_lookup_table[result.query];
	//		const auto homologue = protein_lookup_table[result.target];
	//		assert(cmip != nullptr);
	//		assert(homologue != nullptr);
	//		cmip->addHomologue(mipfinder::Result{homologue, result.bitscore});
	//	}
	//}

	//mipfinder::homology::Results
	//	filterSingleDomainAncestors(const mipfinder::Proteome& proteome,
	//								const mipfinder::homology::Results& results)
	//{
	//	mipfinder::homology::Results filtered;

	//	for (const auto& result : results) {
	//		const auto& ancestor_id = result.target;
	//		const auto ancestor = proteome.find(ancestor_id);
	//		assert(ancestor != nullptr);

	//		unsigned domain_count{0};
	//		for (const auto& interpro_entry : ancestor->interproEntries()) {
	//			if (interpro_entry.type == mipfinder::Interpro::Type::DOMAIN_TYPE ||
	//				interpro_entry.type == mipfinder::Interpro::Type::REPEAT) {
	//				++domain_count;
	//			}
	//		}

	//		/* If it is 0, it is very likely that the annotation of that gene has
	//		* failed since virtually all proteins are made up of small domains.
	//		* If it is exactly 1, we assume that the annotation is correct and it
	//		* really only does have one domain */
	//		if (domain_count != 1) {
	//			filtered.push_back(result);
	//		}
	//	}
	//	return filtered;
	//}

	///* Divides a proteome into potential cMIPs and potential ancestor protein
	// * based on lengths specified in the configuration file. */
	//std::pair<mipfinder::ProteinSet, mipfinder::ProteinSet>
	//	divideProteome(const mipfinder::Proteome& proteome,
	//				   const mipfinder::Configuration& config)
	//{
	//	std::size_t max_mip_length = std::stoi(config["MIP"]["max_mip_length"]);
	//	std::size_t min_ancestor_length = std::stoi(config["MIP"]["min_ancestor_length"]);

	//	const auto cmips = mipfinder::filterByLength(proteome.data(),
	//												 1,
	//												 max_mip_length);
	//	const auto ancestors = mipfinder::filterByLength(proteome.data(),
	//													 min_ancestor_length);

	//	/* Mark respective parts of the proteome as a cMIP or an ancestor */
	//	for (const auto& cmip : cmips) {
	//		assert(cmip->type() == mipfinder::Protein::Type::UNKNOWN);
	//		cmip->setType(mipfinder::Protein::Type::CMIP);
	//	}
	//	for (const auto& ancestor : ancestors) {
	//		assert(ancestor->type() == mipfinder::Protein::Type::UNKNOWN);
	//		ancestor->setType(mipfinder::Protein::Type::ANCESTOR);
	//	}

	//	return std::make_pair(cmips, ancestors);
	//}

	/* Applies a set of filters to results file that represent cMIPS being
	 * searched against potential ancestors. Returns a list of filtered HMMER
	 * results */
	//mipfinder::homology::Results filterAncestors(mipfinder::homology::Results results,
	//										  const mipfinder::Proteome& proteome,
	//										  const mipfinder::Configuration& config)
	//{
	//	/* Apply filters to the unique cMIP vs ancestor results to eliminate
	//	* low-confidence results */
	//	const double ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]);
	//	results = mipfinder::homology::filterByBitscore(results,
	//												 ancestor_bitscore_cutoff);

	//	results = mipfinder::homology::filterByBitscore(results,
	//												 120,
	//												 std::less_equal<>());

	//	/* Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
	//	* protein with a similar domain. This would be a problem for very common
	//	* domains such as kinases, zinc fingers etc. */
	//	const std::size_t hits_to_keep = std::stoi(config["MIP"]["max_ancestor_count"]);
	//	results = mipfinder::homology::keepTopHits(results, hits_to_keep);

	//	//Filter out ancestors that are within 40 a.a of the cMIP
	//	const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
	//	results = mipfinder::homology::filterByLengthDifference(results,
	//														 proteome.data(),
	//														 min_length_diff);

	//	results = filterSingleDomainAncestors(proteome, results);

	//	return results;
	//}

	//void setCmipTypes(const mipfinder::Proteome& proteome)
	//{
	//	for (const auto& protein : proteome.data()) {
	//		if (protein->type() != mipfinder::Protein::Type::CMIP) {
	//			continue;
	//		}
	//		if (protein->homologues().empty()) {
	//			protein->setType(mipfinder::Protein::Type::SINGLE_COPY);
	//		}
	//		else {
	//			protein->setType(mipfinder::Protein::Type::HOMOLOGOUS);
	//		}
	//	}
	//}

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

	//Compare microproteins to each other to establish which ones have homologues (homologous microproteins)
	//and which ones do not (single-copy microproteins).
	template <typename Cont>
	void compareMicroproteinsToMicroproteins(const Cont& potential_microproteins,
											 const mipfinder::Mipfinder::HmmerParameters& parameters,
											 const std::filesystem::path& results_output)
	{
		LOG(INFO) << "Finding homologous relationship between potential microproteins";
		const auto extra_param = "--mx " + parameters.scoring_matrix;
		mipfinder::homology::phmmer(potential_microproteins, potential_microproteins, results_output, extra_param);
		LOG(INFO) << "Finished";
	}


	//Filter out proteins whose existence level hints suggests that they are not translated transcripts
	template <typename Cont>
	requires std::ranges::range<Cont> && requires (Cont c, typename Cont::value_type v)
	{
		typename Cont::value_type;
		&Cont::push_back;
		{ v.existenceLevel() } -> std::convertible_to<std::size_t>;
	}
	Cont removeSpuriousProteins(const Cont& proteome, const std::size_t maximum_allowed_existence_level)
	{
		auto protein_existence_filter = [&](const auto& protein)
		{
			return protein.existenceLevel() <= maximum_allowed_existence_level;
		};
		Cont filtered;
		auto real_proteins = proteome | std::views::filter(protein_existence_filter);
		std::ranges::copy(real_proteins, std::back_inserter(filtered));
		return filtered;
	}

	//Find all potential microproteins from a proteome based on predetermined criteria
	template <typename Cont>
	requires std::ranges::range<Cont> && requires (Cont c, typename Cont::value_type v)
	{
		typename Cont::value_type;
		&Cont::push_back;
		{ v.length() } -> std::convertible_to<std::size_t>;
	}
	Cont findMicroproteins(const Cont& proteome,
						   const mipfinder::Mipfinder::RunParameters& run_params,
						   const mipfinder::Mipfinder::HmmerParameters hmmer_params,
						   const std::filesystem::path& homology_search_results)
	{
		//Filter out all proteins in the proteome that are too long to be microproteins
		auto microprotein_filter = [&](const auto& protein) { return protein.length() <= run_params.maximum_microprotein_length; };
		auto potential_microproteins = proteome | std::views::filter(microprotein_filter);

		//Find homologous microproteins 
		detail::compareMicroproteinsToMicroproteins(potential_microproteins, hmmer_params, homology_search_results);
		auto microprotein_homology_search_results = mipfinder::homology::parseResults(homology_search_results);

		//Filter out all microprotein homology results below @bitscore_cutoff as these do not denote real
		//homologous relationships
		const double lowest_allowed_microprotein_homology_bitscore = hmmer_params.microprotein_homologue_bitscore_cutoff;
		auto true_homologous_microproteins = mipfinder::homology::filterByBitscore(microprotein_homology_search_results, lowest_allowed_microprotein_homology_bitscore);

		//Convert homology search results
		auto true_microproteins = mipfinder::homology::findCorrespondingProteins(true_homologous_microproteins, proteome);


		Cont filtered;
		std::ranges::copy(true_microproteins, std::back_inserter(filtered));
		return filtered;
	}




	template <typename T, typename U>
	requires std::ranges::range<T> && std::ranges::range<U> && std::same_as<std::ranges::range_value_t<U>, mipfinder::Protein>
		std::pair<mipfinder::ProteinList, mipfinder::ProteinList> classifyMicroproteins(T& homology_search_results, U& potential_microproteins)
	{
		/* Create a lookup table of potential microproteins for faster processing */
		std::unordered_map<std::string, mipfinder::Protein*> protein_lookup_table; //Identifiers are keys, Protein objects are values
		//for (auto& potential_microprotein : potential_microproteins) {
		//	protein_lookup_table[potential_microprotein.identifier()] = &potential_microprotein;
		//}

		//Associate potential microproteins with their homologues
		for (auto& result : homology_search_results) {

			////Ignore self-hits, proteins aren't homologous to themselves
			//if (result.query == result.target) {
			//	continue;
			//}

			////If homology search result is not in this group of proteins. This can happen if the
			////homology search results were created using a different list of proteins
			//if (protein_lookup_table.count(result.query) == 0 || protein_lookup_table.count(result.target) == 0) {
			//	continue;
			//}

			//const auto potential_microprotein = std::find(std::begin(protein_lookup_table), , result.query);
			//const auto potential_microprotein_homologue = std::ranges::find(protein_lookup_table, result.target);


			//const auto& potential_microprotein = protein_lookup_table.at(result.query);
			//const auto& potential_microprotein_homologue = protein_lookup_table.at(result.target);

			//if (std::ranges::find(protein_lookup_table, result.query) == std::ranges::end(protein_lookup_table))
			//{
			//	continue;
			//}

			//const mipfinder::homology::Result homologue{.query = potential_microprotein.identifier(), .target = potential_microprotein_homologue.identifier(), .bitscore = result.bitscore};
			//potential_microprotein->addHomologue(mipfinder::Result{.protein = potential_microprotein_homologue, .bitscore = result.bitscore});
		}

		//Find which potential microproteins have homologues (homologous microProteins) and
		//which do not (single-copy)
		mipfinder::ProteinList unique_potential_microproteins;
		mipfinder::ProteinList homologous_potential_microproteins;

		//for (const auto& protein : potential_microproteins) {
		//	if (protein.homologues().empty()) {
		//		unique_potential_microproteins.push_back(protein);
		//	}
		//	else {
		//		homologous_potential_microproteins.push_back(protein);
		//	}
		//}

		return std::make_pair(unique_potential_microproteins, homologous_potential_microproteins);
	}


	//void createFolders()
	//{
	//	//Creates the main results folder for the run
	//	const std::string organism_id = m_run_parameters.organism_identifier;
	//	const std::string folder_name = "results_" + organism_id;

	//	const std::filesystem::path results_folder{folder_name};
	//	std::filesystem::create_directory(results_folder);
	//	results_folder_ = results_folder;

	//	//Create the subfolders for individual software package results
	//	//createResultsFolder() has to be run before this is called, as it sets `results_folder_`
	//	msa_folder_ = results_folder_ / std::filesystem::path{"msa"};
	//	hmmprofile_folder_ = results_folder_ / std::filesystem::path{"hmmprofile"};
	//	homologue_folder_ = results_folder_ / std::filesystem::path{"homologues"};
	//	std::filesystem::create_directory(msa_folder_);
	//	std::filesystem::create_directory(hmmprofile_folder_);
	//	std::filesystem::create_directory(homologue_folder_);
	//}

	std::filesystem::path createMipfinderRunResultsFolder(const std::string organism_identifier)
	{
		//Get current date in the YYYY_MM_DD format
		auto current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		const std::string date_format = "%Y_%m_%d_%H:%M:%S";
		std::stringstream formatted_date;
		formatted_date << std::put_time(std::localtime(&current_time), date_format.c_str());

		//Creates the miPFinder run results folder
		const std::filesystem::path results_folder{formatted_date.str() + "_results_" + organism_identifier};
		std::filesystem::create_directory(results_folder);
		return results_folder;
	}




	template <typename Cont1, typename Cont2>
	void findAncestorHomologuesOfMicroproteins(Cont1& potential_microproteins,
											   Cont2& ancestors,
											   const mipfinder::Mipfinder::HmmerParameters& parameters,
											   const std::filesystem::path& results_output)
	{
		const auto extra_param = "--mx " + parameters.scoring_matrix;
		const std::string extra_phmmer_parameters = "--popen " + std::to_string(parameters.gap_open_probability)
			+ " --pextend " + std::to_string(parameters.gap_extension_probability)
			+ extra_param;

		//mipfinder::proteinToFasta(potential_microproteins, "lol.txt");

		mipfinder::homology::phmmer(potential_microproteins, ancestors, results_output, extra_param);
		//if (config_["DEBUG"]["find_unique_cmip_ancestors"] == "true"
		//	&& config_["DEBUG"]["find_homologous_cmip_ancestors"] == "true") {
		//	
		//}
	}

	//template <typename T>
	//concept HasKeyType = requires
	//{
	//	{ T::key_type } -> std::same_as<typename T::key_type>;
	//};

	//template <typename T>
	//concept HasSubscriptOperator = requires
	//{
	//	{ T::operator[] } -> std::same_as<std::function<typename T::key_type(typename T::key_type)>>;
	//};


	template <typename T>
	requires std::ranges::range<T>
	T keepTopHits(T& t, std::size_t hits_to_keep)
	{
		//static_assert(std::same_as < std::ranges::range_value_t<T>, mipfinder::Protein>);
		//std::unordered_map<mipfinder::Result<mipfinder::Protein, mipfinder::Protein>, std::string> count_table;
		//std::unordered_map<mipfinder::Result<int, int>, std::string> m;
		//std::unordered_map<std::ranges::range_value_t<T>, std::size_t> count_table;
		T filtered;
		//for (const auto& elem : t) {
		//	count_table[elem] += 1;
		//	if (count_table[elem] <= hits_to_keep) {
		//		filtered.push_back(elem);
		//	}
		//}
		return filtered;
	}

}



mipfinder::Mipfinder::Mipfinder(const std::filesystem::path& configuration_file)
{
	configureLogger();

	//Set up configuration parameters from file
	LOG(DEBUG) << "Setting up miPFinder run configuration";
	//extractConfiguration(configuration_file);
	Configuration config{configuration_file};
	m_hmmer_parameters = HmmerParameters{
		.microprotein_homologue_bitscore_cutoff = std::stod(config["HMMER"]["matrix"]),
		.ancestor_bitscore_cutoff = std::stod(config["HMMER"]["ancestor_bitscore_cutoff"]),
		.gap_open_probability = std::stod(config["HMMER"]["gap_open_probability"]),
		.gap_extension_probability = std::stod(config["HMMER"]["gap_extend_probability"]),
		.scoring_matrix = config["HMMER"]["matrix"]
	};

	m_file_parameters = FileParamaters{
		.input_proteome = config["TARGET"]["input_proteome"],
		.known_microprotein_list = config["TARGET"]["known_micropotein_list"],
		.interpro_database = config["INTERPRO"]["interpro_database"],
		.gene_ontology_database = config["GO"]["go_database"],
		.uniprot_to_intepro_id_conversion_file = config["INTERPRO"]["uniprot_to_interpro_id_conversion"],
		.uniprot_to_go_id_conversion_file = config["GO"]["uniprot_to_go_id_conversion"],
	};

	m_run_parameters = RunParameters{
		.maximum_microprotein_length = std::stoul(config["MIP"]["maximum_microprotein_length"]),
		.minimum_ancestor_length = std::stoul(config["MIP"]["minimum_ancestor_length"]),
		.maximum_homologues_per_microprotein = std::stoul(config["MIP"]["maximum_homologues_per_microprotein"]),
		.minimum_length_difference = std::stoul(config["MIP"]["minimum_length_difference"]),
		.maximum_ancestor_count = std::stoul(config["MIP"]["maximum_ancestor_count"]),
		.maximum_protein_existence_level = std::stoul(config["TARGET"]["maximum_protein_existence_level"]),
		.output_format = config["REPORT"]["format"],
		.organism_identifier = config["TARGET"]["organism_identifier"],
	};

	LOG(DEBUG) << "Checking file dependencies";
	checkFileDependencies(m_file_parameters);

	//LOG(DEBUG) << "Setting up GO database";
	//mipfinder::Go go_database{config_["GO"]["go_database"]};
	//LOG(DEBUG) << "Adding GO ids to proteins";
	//addGoIdentifiers(proteome_,
	//				 go_database,
	//				 config_["GO"]["uniprot_to_go"]);

	//LOG(DEBUG) << "Setting up InterPro database";
	//mipfinder::Interpro interpro_database{config_["INTERPRO"]["interpro_database"]};
	//LOG(DEBUG) << "Adding InterPro ids to proteins";
	//addInterproIdentifiers(proteome_,
	//					   interpro_database,
	//					   config_["INTERPRO"]["uniprot_to_interpro"]);

	LOG(DEBUG) << "Creating miPFinder run folders";
	m_results_folder = detail::createMipfinderRunResultsFolder(m_run_parameters.organism_identifier);


	////WIP: DO this at the end
	///* Provide a copy of the run parameters with the results */
	//std::filesystem::path config_file_copy = results_folder_ / configuration_file;
	//std::filesystem::copy(configuration_file, config_file_copy, std::filesystem::copy_options::overwrite_existing);
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
		LOG(INFO) << "Loading proteome...";
		auto proteome = detail::loadProteome(m_file_parameters.input_proteome);
		LOG(INFO) << "Detected " << proteome.size() << " proteins";

		//Remove proteins whose existence level implies their transcripts do not
		const std::size_t maximum_allowed_existence_level = m_run_parameters.maximum_protein_existence_level;
		auto real_proteins = detail::removeSpuriousProteins(proteome, maximum_allowed_existence_level);

		//Find all potential microproteins in the proteome based on size
		//-----------------------------
		const std::filesystem::path classified_microproteins = m_results_folder / "all_cmips_vs_cmips.txt";
		auto potential_microproteins = detail::findMicroproteins(proteome, m_run_parameters, m_hmmer_parameters, classified_microproteins);



		//Find all potential ancestors in the proteome based on size
		//const potential_ancestors = detail::findAncestors(proteome);

		const std::size_t minimum_allowed_ancestor_length = m_run_parameters.minimum_ancestor_length;
		auto ancestor_filter = [&](const auto& protein) { return protein.length() <= minimum_allowed_ancestor_length; };
		auto potential_ancestors = real_proteins | std::views::filter(ancestor_filter);





















		///* Deal with single copy cMIPS */
		///* Take all single-copy cMIPS and phmmer them against large proteins to find
		// * their ancestors */
		//auto [unique_potential_microproteins, homologous_unique_microproteins] = detail::classifyMicroproteins(high_confidence_microprotein_homologues, potential_microproteins);

		////Find all potential ancestors in the proteome based on size


		//auto unique_microproteins_vs_ancestors = m_results_folder / "unique_vs_ancestor.txt";
		//detail::findAncestorHomologuesOfMicroproteins(unique_potential_microproteins, potential_ancestors, m_hmmer_parameters, unique_microproteins_vs_ancestors);
		//auto unique_microprotein_vs_ancestor_homology_search_results = mipfinder::homology::parseResults(unique_microproteins_vs_ancestors);

		///* Apply filters to the unique cMIP vs ancestor results to eliminate
		// * low-confidence results */
		//const auto maximum_allowed_ancestor_homology_bitscore = m_hmmer_parameters.ancestor_bitscore_cutoff;
		//const auto minimum_allowed_ancestor_homology_bitscore = 120;
		//auto ancestor_bitscore_filter = [&](const auto& homology_search_result)
		//{
		//	return homology_search_result.bitscore <= maximum_allowed_ancestor_homology_bitscore || homology_search_result.bitscore >= minimum_allowed_ancestor_homology_bitscore;
		//};
		//auto high_confidence_ancestors = unique_microprotein_vs_ancestor_homology_search_results | std::views::filter(ancestor_bitscore_filter);

		///* Keep the top x ancestors for each MIP to ensure we don't pick up EVERY
		//* protein with a similar domain. This would be a problem for very common
		//* domains such as kinases, zinc fingers etc. */
		//const auto maximum_homologous_ancestors_per_microprotein_to_keep = m_run_parameters.maximum_homologues_per_microprotein;
		//auto filtered_ancestor_homologues = detail::keepTopHits(high_confidence_ancestors, maximum_homologous_ancestors_per_microprotein_to_keep);





//		const std::size_t hits_to_keep = std::stoi(config["MIP"]["max_ancestor_count"]);
//		results = mipfinder::homology::keepTopHits(results, hits_to_keep);
//
//		//Filter out ancestors that are within 40 a.a of the cMIP
//		const int min_length_diff = std::stoi(config["MIP"]["min_length_difference"]);
//		results = mipfinder::homology::filterByLengthDifference(results,
//															 proteome.data(),
//															 min_length_diff);
//
//		results = filterSingleDomainAncestors(proteome, results);
//
//		return results;



		//const auto unique_vs_ancestor_file =
		//	phmmerCmipsVsAncestors(unique_cmips, "unique_vs_proteome.txt");

		//auto unique_vs_ancestor_results = hmmer::parseResults(unique_vs_ancestor_file);

		//unique_vs_ancestor_results = filterAncestors(unique_vs_ancestor_results,
		//											 proteome_,
		//											 config_);

		//associateAncestorsWithCmips(proteome_,
		//							unique_vs_ancestor_results);

		///* Deal with homologous cMIPS */
		///* Build HMM profiles of all homologous cMIPS and hmmsearch these profiles
		// * against the ancestors */
		//writeHomologuesToFasta(homologous_cmips);
		//alignHomologousCmips(homologue_folder_);

		//if (config_["DEBUG"]["create_hmm_profiles"] == "true") {
		//	hmmer::createHmmProfiles(msa_folder_, hmmprofile_folder_);
		//}

		//std::filesystem::path filename{"all_hmmprofiles.txt"};
		//std::filesystem::path unified_hmmprofile_file = results_folder_ / filename;
		//const auto all_hmmer_profiles = hmmer::concatenateHmmProfiles(hmmprofile_folder_,
		//															  unified_hmmprofile_file);

		///* Merges all individual hmm profile files into one large one */
		//const auto homologous_vs_ancestor_search =
		//	hmmsearchHomologousCmips(all_hmmer_profiles);

		//auto homologous_vs_ancestor_results = hmmer::parseResults(homologous_vs_ancestor_search);

		//homologous_vs_ancestor_results = filterAncestors(homologous_vs_ancestor_results,
		//												 proteome_,
		//												 config_);

		//associateAncestorsWithCmips(proteome_,
		//							homologous_vs_ancestor_results);

		//setCmipTypes(proteome_);

		///* Assign KNOWN_MIP type to all known mips */
		//Proteome known_mips(config_["MIP"]["known_mips_fasta"]);
		//for (const auto& mip : known_mips.data()) {
		//	auto mip_in_proteome = proteome_.find(mip->identifier());
		//	if (mip_in_proteome) {
		//		mip_in_proteome->setType(mipfinder::Protein::Type::KNOWN_MIP);
		//	}
		//}

		///* Filter out proteins with more than @maximum_homologues_allowed homologues */
		//mipfinder::ProteinSet proteins_to_score;
		//for (const auto& protein : proteome_.data()) {
		//	const double max_homologues_allowed = std::stod(config_["MIP"]["maximum_homologues"]);
		//	if (protein->homologues().size() > max_homologues_allowed) {
		//		continue;
		//	}
		//	proteins_to_score.insert(protein);
		//}

		//assignHmmerScores(proteins_to_score, &scoring_algorithm);
		//assignHmmerScores(proteins_to_score, &scoring_algorithm);

		///* Filter out proteins with a 0 score */
		//mipfinder::ProteinSet non_zero;
		//for (const auto& protein : proteome_.data()) {
		//	if (protein->score() != 0) {
		//		non_zero.insert(protein);
		//	}
		//}

		///* Filter out proteins that are not CMIPS */
		//mipfinder::ProteinSet only_cmip_type;
		//for (const auto& protein : non_zero) {
		//	if ((protein->type() == mipfinder::Protein::Type::SINGLE_COPY) ||
		//		(protein->type() == mipfinder::Protein::Type::HOMOLOGOUS) ||
		//		(protein->type() == mipfinder::Protein::Type::KNOWN_MIP)) {
		//		only_cmip_type.insert(protein);
		//	}
		//}

		///* Sort proteome based on score */
		//auto comp = [](const mipfinder::Protein* lhs, const mipfinder::Protein* rhs)
		//{
		//	return lhs->score() > rhs->score();
		//};
		//std::vector<Protein*> sorted_proteome{only_cmip_type.begin(), only_cmip_type.end()};
		//std::sort(sorted_proteome.begin(), sorted_proteome.end(), comp);

		//std::filesystem::path results_fasta_file{"cmip_results.fasta"};
		//std::filesystem::path results_fasta_location = results_folder_ / results_fasta_file;
		//mipfinder::proteinToFasta(sorted_proteome, results_fasta_location);

		//LOG(INFO) << "Writing final report";
		//const std::string report_format = config_["REPORT"]["format"];
		//std::filesystem::path final_results = results_folder_ / "final_results.txt";
		//mipfinder::printer::createReport(report_format,
		//								 '\t',
		//								 sorted_proteome,
		//								 final_results);

		//LOG(INFO) << "mpf v2.0 has finished.";
	}

	//void Mipfinder::createFolders()
	//{
	//	//Creates the main results folder for the run
	//	const std::string organism_id = m_run_parameters.organism_identifier;
	//	const std::string folder_name = "results_" + organism_id;

	//	const std::filesystem::path results_folder{folder_name};
	//	std::filesystem::create_directory(results_folder);
	//	results_folder_ = results_folder;

	//	//Create the subfolders for individual software package results
	//	//createResultsFolder() has to be run before this is called, as it sets `results_folder_`
	//	msa_folder_ = results_folder_ / std::filesystem::path{"msa"};
	//	hmmprofile_folder_ = results_folder_ / std::filesystem::path{"hmmprofile"};
	//	homologue_folder_ = results_folder_ / std::filesystem::path{"homologues"};
	//	std::filesystem::create_directory(msa_folder_);
	//	std::filesystem::create_directory(hmmprofile_folder_);
	//	std::filesystem::create_directory(homologue_folder_);
	//}

	//std::filesystem::path Mipfinder::phmmerAgainstSelf(const mipfinder::ProteinSet& cmips)
	//{
	//	std::filesystem::path results_file{"all_cmips_against_cmips.txt"};
	//	const auto hmmer_cmips_vs_cmips = results_folder_ / results_file;

	//	if (config_["DEBUG"]["phmmer_cmips"] == "false") {
	//		return std::filesystem::path{hmmer_cmips_vs_cmips};
	//	}

	//	LOG(INFO) << "Grouping cMIPs based on homology";

	//	const auto scoring_matrix = config_["HMMER"]["matrix"];
	//	const auto extra_param = " --mx " + scoring_matrix;

	//	mipfinder::homology::phmmer(cmips, cmips, hmmer_cmips_vs_cmips.string(), extra_param);
	//	LOG(INFO) << "Finished grouping cMIPS";
	//	return hmmer_cmips_vs_cmips;
	//}

	//void Mipfinder::writeHomologuesToFasta(const mipfinder::ProteinSet& homologous_cmips)
	//{
	//	//if (config_["DEBUG"]["write_homologous_groups_fasta"] == "false") { return; }

	//	LOG(INFO) << "Creating FASTA files of homologous cMIPs";
	//	int grouped_fasta_files{0};
	//	for (const auto& protein : homologous_cmips) {
	//		mipfinder::ProteinSet grouped_homologues;
	//		grouped_homologues.insert(protein);

	//		for (const auto& homologue : protein->homologues()) {
	//			grouped_homologues.insert(homologue.protein);
	//		}

	//		const auto id = protein->identifier();

	//		const std::filesystem::path unaligned_filename{id + "_homologues.fasta"};
	//		const std::filesystem::path unaligned_file_location{homologue_folder_
	//			/ unaligned_filename};

	//		proteinToFasta(grouped_homologues, unaligned_file_location.string());
	//		++grouped_fasta_files;
	//	}
	//	LOG(DEBUG) << "Created " << grouped_fasta_files << " homologous cMIP sequence "
	//		<< "FASTA files";
	//}

	//void Mipfinder::alignHomologousCmips(const std::filesystem::path& unaligned_seq_folder)
	//{
	//	//if (config_["DEBUG"]["align_homologous_cmips"] == "false") { return; }

	//	LOG(INFO) << "Aligning groups of homologous cMIPs";

	//	typedef std::filesystem::directory_iterator DirectoryIter;

	//	int msas_created{0};
	//	for (const auto& file : DirectoryIter{unaligned_seq_folder}) {

	//		//Unaligned FASTA files are named as `PROTEIN_ID`_homologues.fasta. We
	//		//need to extract the `PROTEIN_ID` part for the aligned filename.
	//		const auto unaligned_filename = file.path().filename().string();
	//		const auto tokens = mipfinder::tokenise(unaligned_filename, '_');
	//		const auto protein_id = tokens[0];

	//		const std::filesystem::path unaligned_file_path = file.path();

	//		const std::filesystem::path msa_filename{protein_id + "_aligned.fasta"};
	//		const std::filesystem::path msa_file_full_path =
	//			msa_folder_ / msa_filename;

	//		const std::string clustalo_command = "clustalo -i "
	//			+ unaligned_file_path.string()
	//			+ " -o "
	//			+ msa_file_full_path.string();

	//		int sys_call_result = std::system(clustalo_command.c_str());
	//		if (sys_call_result != 0) {
	//			continue;
	//		}
	//		++msas_created;
	//	}
	//	LOG(INFO) << "Created " << msas_created << " MSA from grouped cMIP homologues";
	//}

	//std::filesystem::path
	//	Mipfinder::hmmsearchHomologousCmips(const std::filesystem::path& hmmprofiles_file)
	//{
	//	unsigned int min_ancestor_length =
	//		std::stoi(config_["MIP"]["min_ancestor_length"]);

	//	const auto ancestors =
	//		mipfinder::filterByLength(proteome_.data(), min_ancestor_length);

	//	const std::filesystem::path results_filename{"homologous_vs_proteome.txt"};
	//	const std::filesystem::path results_file_location =
	//		results_folder_ / results_filename;

	//	if (config_["DEBUG"]["find_homologous_cmip_ancestors"] == "true") {
	//		mipfinder::homology::hmmsearch(hmmprofiles_file, ancestors, results_file_location);
	//	}
	//	return results_file_location;
	//}

	//void Mipfinder::assignHmmerScores(const mipfinder::ProteinSet& proteins,
	//								  HmmerScoringAlgorithm algorithm)
	//{
	//	int score_counter{0};
	//	for (const auto& protein : proteins) {
	//		/* Only cMIPs have ancestors */
	//		if (protein->ancestors().empty()) {
	//			continue;
	//		}

	//		for (const auto& ancestor : protein->ancestors()) {
	//			const double bitscore = ancestor.bitscore;

	//			double score = algorithm(HmmerScoringData{protein, ancestor.protein, bitscore});
	//			protein->changeScore(score);
	//			++score_counter;
	//		}
	//	}
	//}

}
