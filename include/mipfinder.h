#ifndef MIPFINDER_MIPFINDER_H
#define MIPFINDER_MIPFINDER_H

#include <filesystem>
#include <string>

#include "aliases.h"
#include "configuration.h"
#include "hmmer_scoring_data.h"
#include "proteome.h"

namespace mipfinder {
  class Mipfinder {
  public:
    Mipfinder(const std::filesystem::path& configuration_file);
    
    //It is responsible for running the actual mipfinder pipeline. See 
    //documentation for details. 
    void run();

    void writeOutput(std::string filename);

  private:  
    Configuration config_;
    Proteome proteome_;

    /* All folders that mipfinder needs to run properly */
    
    std::filesystem::path results_folder_;
    std::filesystem::path msa_folder_;
    std::filesystem::path hmmprofile_folder_;
    std::filesystem::path homologue_folder_;

    //Creates all necessary results folder to run mipfinder v2.0
    void createFolders();

    //Uses phmmer to compare each cMIP against each other to find cMIPs with no
    //identified homologues (unique cMIPS) and groups of homologous cMIPS.
    //Returns a path to the phmmer output file.
    std::filesystem::path phmmerAgainstSelf(const mipfinder::ProteinSet& cmips);

    std::filesystem::path
    phmmerCmipsVsAncestors(const mipfinder::ProteinSet& cmips,
                           const std::filesystem::path& results_filename);

    //Creates a FASTA file of all homologous cMIP sequences to be used as an 
    //input for clustalo. For each protein it creates a file called 
    //`PROTEIN_ID`_homologues.fasta in `ORGANISM_NAME`/msa/ folder, where
    //`ORGANISM_NAME` is specified in the configuration and `PROTEIN_ID` is the
    //protein UniProt identifier.
    void writeHomologuesToFasta(const mipfinder::ProteinSet& homologous_cmips);

    //Runs clustalo on each file in the specified directory.
    void alignHomologousCmips(const std::filesystem::path& unaligned_seq_dir);

    //Performs a hmmsearch on the proteome using a file containing one or more
    //hmmprofiles. Returns a path to the hmmsearch results file which is in a 
    //table format
    std::filesystem::path
    hmmsearchHomologousCmips(const std::filesystem::path& hmmprofiles_file);

    typedef double (*HmmerScoringAlgorithm)(const HmmerScoringData&);

    void assignHmmerScores(const mipfinder::ProteinSet& proteome,
                           HmmerScoringAlgorithm algorithm);
  };
}

#endif
