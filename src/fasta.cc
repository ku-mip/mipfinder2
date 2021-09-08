#include <iostream>
#include <filesystem>

#include "fasta.h"

mipfinder::fasta::Entries mipfinder::fasta::parse(const std::filesystem::path& file)
{
    std::ifstream stream;
    stream.open(file);
    if (stream.is_open()) {
        return mipfinder::fasta::parse(stream);
    }
    else {
        throw std::runtime_error("Could not open " + file.string() + ", aborting");
    }
}

mipfinder::fasta::Entries mipfinder::fasta::parse(std::ifstream& stream)
{
    /* This ensures we are only parsing lines that are part of a FASTA record */
    bool fasta_record_found = false;

    std::string fasta_header;
    std::string fasta_sequence;

    mipfinder::fasta::Entries all_fasta_records;

    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty()) {
            continue;
        }

        //If the FASTA header contains no information other than the '>' character marking it as a header,
        //ignore the line.
        if (line.front() == '>' && line.length() == 1) {
            continue;
        }

        if ((line.front() == '>') && (!fasta_record_found)) {
            fasta_record_found = true;
            fasta_header = line;
        }

        else if (!(line.front() == '>') && fasta_record_found) {
            fasta_sequence += line;
        }

        else if ((line.front() == '>') && fasta_record_found) {
            /* If we have found a new record but the current record did not have a
             * sequence, ignore it */
            if (fasta_sequence.empty()) {
                fasta_sequence.clear();
                fasta_header = line;
                continue;
            }
            else {
                all_fasta_records.emplace_back(mipfinder::fasta::Entry{ .header = fasta_header, .sequence = fasta_sequence });

                /* Reset the variable contents */
                fasta_sequence.clear();
                fasta_header = line;
            }
        }

        else if (!fasta_record_found) {
            /* If no record is found or being processed, the data is bad */
            continue;
        }
    }

    /* If final sequence is empty, ignore the header */
    if (fasta_sequence.empty()) {
        return all_fasta_records;
    }

    all_fasta_records.emplace_back(mipfinder::fasta::Entry{ .header = fasta_header, .sequence = fasta_sequence });
    return all_fasta_records;
}
