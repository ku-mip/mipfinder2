#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

#include "configuration.h"
#include "easylogging++.h"
#include "file.h"

#include "hmmer.h"
#include "interpro.h"
#include "protein.h"
#include "proteome.h"

namespace
{
	mipfinder::Interpro::Entries
		filterInterproByType(const mipfinder::Interpro::Entries& entries,
							 const mipfinder::Interpro::Type& type)
	{
		mipfinder::Interpro::Entries filtered;
		for (const auto& entry : entries) {
			if (entry.type == type) {
				filtered.push_back(entry);
			}
		}
		return filtered;
	}

}


namespace mipfinder::hmmer
{
	template <>
	void phmmer(const std::filesystem::path& query_file,
				const std::filesystem::path& database_file,
				const std::filesystem::path& results_file,
				const std::string& extra_parameters)
	{
		std::string phmmer_output = "hmmer_output.fasta";
		std::string phmmer_output_table = "hmmer_output_table.fasta";

		std::vector<std::string> command_tokens{std::string{"phmmer"},
												std::string{"-o /dev/null"},
												std::string{"--tblout " + results_file.string()},
												extra_parameters,
												query_file.string(),
												database_file.string()};

		std::string phmmer_command;
		for (const auto& token : command_tokens) {
			const char token_separator{' '};
			phmmer_command += token + token_separator;
		}

		LOG(DEBUG) << "Starting phmmer with the following command: \"" << phmmer_command << "\"";
		int sys_call_result = std::system(phmmer_command.c_str());
		if (sys_call_result != 0) {
			LOG(ERROR) << "Failed to call a subprocess in phmmer()";
		}
	}

	void buildHmmerProfile(const std::filesystem::path& msa_file,
						   const std::filesystem::path& output_file,
						   const std::string& extra_parameters)
	{
		const std::vector<std::string> command_tokens{std::string{"hmmbuild"},
													  std::string{"-o /dev/null"},
													  extra_parameters,
													  output_file.string(),
													  msa_file.string()};

		std::string hmmbuild_command;
		for (const auto token : command_tokens) {
			hmmbuild_command += token + " "; //Space to separate the tokens
		}

		int sys_call_result = std::system(hmmbuild_command.c_str()); //std::system expects a C-style string
		if (sys_call_result != 0) {
			LOG(ERROR) << "Failed to call a subprocess in builHmmerProfile()";
		}
	}

	void hmmsearch(const std::filesystem::path& profile_file,
				   const mipfinder::ProteinList& database,
				   const std::filesystem::path& results_file,
				   const std::string& extra_parameters)
	{
		const std::filesystem::path results_path = results_file.parent_path();

		std::filesystem::path database_file{"hmmsearch_database.txt"};
		std::filesystem::path database_file_location = results_path / database_file;

		mipfinder::proteinToFasta(database, database_file_location);

		hmmsearch(profile_file,
				  database_file_location,
				  results_file,
				  extra_parameters);
	}

	void hmmsearch(const std::filesystem::path& profile_file,
				   const std::filesystem::path& database_file,
				   const std::filesystem::path& results_file,
				   const std::string& extra_parameters)
	{
		const std::vector<std::string> command_tokens{std::string{"hmmsearch"},
													  std::string{"-o /dev/null"},
													  std::string{"--tblout"},
													  results_file.string(),
													  profile_file.string(),
													  database_file.string(),
													  extra_parameters};

		std::string hmmsearch_command;
		for (const auto token : command_tokens) {
			hmmsearch_command += token + " "; //Space to separate the tokens
		}

		LOG(INFO) << "Starting hmmsearch with the following command: \"" << hmmsearch_command << "\"";
		int sys_call_result = std::system(hmmsearch_command.c_str()); //std::system expects a C-style string
		if (sys_call_result != 0) {
			LOG(ERROR) << "Failed to call a subprocess in hmmsearch()";
		}
	}

	mipfinder::hmmer::Results parseResults(const std::filesystem::path& results_file)
	{
		auto stream = mipfinder::file::open(results_file);

		mipfinder::hmmer::Results results;
		std::string line;
		while (std::getline(stream, line)) {
			if (line.front() == '#') { //'#' lines are comments
				continue;
			}

			const auto tokens = mipfinder::tokenise(line, ' ');
			const std::string query = tokens[2];
			const std::string target = tokens[0];
			double bitscore = stod(tokens[5]);

			results.emplace_back(mipfinder::hmmer::Result{.query = query, .target = target, .bitscore = bitscore});
		}
		return results;
	}

	mipfinder::hmmer::Results
		keepTopHits(const mipfinder::hmmer::Results& results, std::size_t hits_to_keep)
	{
		std::unordered_map<std::string, unsigned int> count_table;
		mipfinder::hmmer::Results filtered;
		for (const auto& result : results) {
			count_table[result.query] += 1;
			if (count_table[result.query] <= hits_to_keep) {
				filtered.push_back(result);
			}
		}
		return filtered;
	}


	mipfinder::hmmer::Results
		removeSelfHits(const mipfinder::hmmer::Results& results)
	{
		mipfinder::hmmer::Results filtered;
		for (const auto& result : results) {
			if (result.query == result.target) {
				continue;
			}
			filtered.push_back(result);
		}
		return filtered;
	}

	void createHmmProfiles(const std::filesystem::path& msa_dir,
						   const std::filesystem::path& output_dir)
	{
		typedef std::filesystem::directory_iterator DirectoryIter;

		int profiles_built{0};
		for (const auto& msa_file : DirectoryIter{msa_dir}) {
			//MSA files are named as `PROTEIN_ID`_aligned.fasta. We
			//need to extract the `PROTEIN_ID` part for the aligned filename.
			const auto msa_filename = msa_file.path().filename().string();
			const auto tokens = mipfinder::tokenise(msa_filename, '_');
			const auto protein_id = tokens[0];

			const auto msa_file_path = msa_file.path();

			const std::filesystem::path hmmprofile_filename{protein_id + "_profile.hmm"};
			const std::filesystem::path hmmprofile_file_location =
				output_dir / hmmprofile_filename;

			//Ensures that the HMMER profile name matches protein_id, rather than using
			//the filename.
			const auto hmmprofile_name_command = "-n " + protein_id;
			mipfinder::hmmer::buildHmmerProfile(msa_file_path.string(),
												hmmprofile_file_location,
												hmmprofile_name_command);
			++profiles_built;
		}

	}

	std::filesystem::path
		concatenateHmmProfiles(const std::filesystem::path& hmmer_profile_dir,
							   const std::filesystem::path& output_file)
	{
		typedef std::filesystem::directory_iterator DirectoryIter;

		std::ofstream of;
		of.open(output_file, std::ios::trunc);

		for (const auto& profile_file : DirectoryIter{hmmer_profile_dir}) {
			std::ifstream f;
			f.open(profile_file.path());
			std::string line;
			while (getline(f, line)) {
				of << line << "\n";
			}
		}
		return output_file;
	}

	mipfinder::hmmer::Results
		filterByLengthDifference(const mipfinder::hmmer::Results& results,
								 const mipfinder::ProteinSet& proteins,
								 unsigned int min_length_difference)
	{
		mipfinder::hmmer::Results filtered;

		std::unordered_map<std::string, mipfinder::Protein*> lookup_table;
		for (const auto& protein : proteins) {
			lookup_table[protein->identifier()] = protein;
		}

		for (const auto result : results) {
			const auto query = lookup_table.at(result.query);
			const auto target = lookup_table.at(result.target);

			if (target->length() <= query->length()) {
				continue;
			}
			assert(target->length() > query->length());

			const auto length_difference = target->length() - query->length();
			if (length_difference < min_length_difference) {
				continue;
			}

			filtered.push_back(result);
		}
		return filtered;
	}

	mipfinder::hmmer::Results
		filterByProteinFamilyAndDomain(const mipfinder::hmmer::Results& results,
									   const mipfinder::Proteome& proteome)
	{
		mipfinder::hmmer::Results filtered_results;
		for (const auto& result : results) {
			const auto query = proteome.find(result.query);
			const auto target = proteome.find(result.target);

			const auto query_interpro_entries = query->interproEntries();
			const auto target_interpro_entries = target->interproEntries();

			const auto query_family = filterInterproByType(query_interpro_entries,
														   mipfinder::Interpro::Type::FAMILY);

			const auto query_domain = filterInterproByType(query_interpro_entries,
														   mipfinder::Interpro::Type::DOMAIN_TYPE);
			mipfinder::Interpro::Entries query_merged;
			for (const auto& i : query_family) {
				query_merged.push_back(i);
			}
			for (const auto& i : query_domain) {
				query_merged.push_back(i);
			}

			const auto target_family = filterInterproByType(target_interpro_entries,
															mipfinder::Interpro::Type::FAMILY);

			const auto target_domain = filterInterproByType(target_interpro_entries,
															mipfinder::Interpro::Type::DOMAIN_TYPE);
			mipfinder::Interpro::Entries target_merged;
			for (const auto& i : target_family) {
				target_merged.push_back(i);
			}
			for (const auto& i : target_domain) {
				target_merged.push_back(i);
			}

			for (const auto& interpro_entry : query_merged) {
				if (std::find(target_merged.begin(),
					target_merged.end(),
					interpro_entry) == query_merged.end()) {
					continue;
				}
				filtered_results.push_back(result);
				break;
			}
		}
		return filtered_results;
	}

	mipfinder::hmmer::Results
		filterByHomologueCount(const mipfinder::hmmer::Results& results,
							   unsigned int maximum_homologues)
	{
		/* Map proteins to their identified homologues */
		std::unordered_map<std::string, std::vector<std::string>> lookup_table;
		for (const auto& result : results) {
			lookup_table[result.query].push_back(result.target);
		}

		mipfinder::hmmer::Results filtered_results;
		for (const auto& result : results) {
			if (lookup_table[result.query].size() > maximum_homologues) {
				continue;
			}
			filtered_results.push_back(result);
		}
		return filtered_results;
	}

	mipfinder::ProteinSet
		convertToProtein(const mipfinder::hmmer::Results& results,
						 const mipfinder::ProteinSet& proteins)
	{
		/* Lookup table for fast searching */
		std::unordered_map<std::string, mipfinder::Protein*> lookup_table;
		for (const auto& protein : proteins) {
			lookup_table[protein->identifier()] = protein;
		}

		mipfinder::ProteinSet found_proteins;
		for (const auto& result : results) {
			const auto query_protein = lookup_table.at(result.query);
			const auto ancestor_protein = lookup_table.at(result.target);

			assert(query_protein != nullptr);
			assert(ancestor_protein != nullptr);

			found_proteins.insert(ancestor_protein);
		}

		return found_proteins;
	}


}
