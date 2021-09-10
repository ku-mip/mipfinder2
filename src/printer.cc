#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <iostream>
#include <fstream>

#include "helpers.h"
#include "printer.h"
#include "protein.h"

namespace
{
  std::string createHeader(const std::string& format, char delimiter)
  {
    const auto format_tokens = mipfinder::split(format, '%');
    std::string header;

    for (const auto& token : format_tokens) {
      if (token == format_tokens.back()) {
        delimiter = ' ';
      }

      if (token == R"(%acc)") {
        header = header + "ACCESSION_NAME" + delimiter;
      }
      else if (token == "%length") {
        header = header + "LENGTH" + delimiter;
      }
      else if (token == "%ancestors") {
        header = header + "ANCESTORS" + delimiter;
      }
      else if (token == "%score") {
        header = header + "SCORE" + delimiter;
      }
      else if (token == "%homologues") {
        header = header + "HOMOLOGUES" + delimiter;
      }
      else if (token == "%description") {
        header = header + "DESCRIPTION" + delimiter;
      }
      else if (token == "%gene_ontology") {
        header = header + "GENE_ONTOLOGY" + delimiter;
      }
      else if (token == "%type") {
        header = header + "MIP_TYPE" + delimiter;
      }
      else if (token == "%ancestor_nr") {
        header = header + "NR_OF_ANCESTORS" + delimiter;
      }
      else if (token == "%interpro") {
        header = header + "INTERPRO_DOMAINS" + delimiter;
      }
      else if (token == "%ancestor_domains") {
        header = header + "ANCESTOR_INTERPRO_DOMAINS" + delimiter;
      }
      else if (token == "%ancestor_go") {
        header = header + "ANCESTOR_GO_DOMAINS" + delimiter;
      }
      else if (token == "%instability") {
        header = header + "INSTABILITY_INDEX" + delimiter;
      }
      else if (token == "%combined_score") {
        header = header + "INSTABILITY_INDEX*SCORE" + delimiter;
      }
      else {
        continue;
      }
    }
    return header;
  }

}
//
//namespace mipfinder::printer
//{
//  void createReport(const std::string& format,
//                    char delimiter,
//                    const std::vector<Protein*>& proteins,
//                    const std::filesystem::path& output_file)
//  {
//
//    std::ofstream of;
//    of.open(output_file);
//    const auto header = createHeader(format, delimiter);
//    of << header << "\n";
//
//    for (const auto& protein : proteins) {
//      std::string report_line;
//      auto format_tokens = mipfinder::split(format, '%');
//      for (const auto& token : format_tokens) {
//        if (token == "%acc") {
//          report_line = report_line + protein->identifier().to_string() + delimiter;
//        }
//
//        else if (token == "%length") {
//          report_line = report_line + std::to_string(protein->length()) + delimiter;
//        }
//
//        //else if (token == "%ancestors") {
//        //  std::string all_ancestor_names;
//        //  for (const auto& ancestor : protein->ancestors()) {
//        //    const auto bitscore = ancestor.bitscore;
//        //    std::stringstream stream;
//        //    stream << std::fixed << std::setprecision(1) << bitscore;
//        //    const auto bitscore_string = stream.str();
//
//        //    all_ancestor_names = all_ancestor_names + ancestor.protein->identifier()
//        //                         + "(" + bitscore_string + ");";
//        //  }
//        //  report_line = report_line + all_ancestor_names + delimiter;
//        //}
//
//        //else if (token == "%score") {
//        //  report_line = report_line + std::to_string(protein->score()) + delimiter;
//        //}
//
//        //else if (token == "%homologues") {
//        //  std::string homologue_names;
//        //  for (const auto& homologue_entry : protein->homologues()) {
//        //    const auto bitscore = homologue_entry.bitscore;
//        //    std::stringstream stream;
//        //    stream << std::fixed << std::setprecision(1) << bitscore;
//        //    const auto bitscore_string = stream.str();
//
//        //    homologue_names = homologue_names + homologue_entry.protein->identifier()
//        //                      + "(" + bitscore_string + ");";
//        //  }
//        //  report_line = report_line + homologue_names + delimiter;
//        //}
//
//        else if (token == "%description") {
//          report_line = report_line + protein->description() + delimiter;
//        }
//
//        //else if (token == "%gene_ontology") {
//        //  std::string go_terms;
//        //  for (const auto& record : protein->goEntries()) {
//        //    go_terms = go_terms + record.process_name + ";";
//        //  }
//        //  report_line = report_line + go_terms + delimiter;
//        //}
//
//        //else if (token == "%type") {
//        //  report_line = report_line + protein->type_to_string() + delimiter;
//        //}
//
//        //else if (token == "%instability") {
//        //  report_line = report_line + std::to_string(mipfinder::instability_index(protein->sequence())) + delimiter;
//        //}
//        //else if (token == "%ancestor_nr") {
//        //  report_line = report_line + std::to_string(protein->ancestors().size()) + delimiter;
//        //}
//        //else if (token == "%interpro") {
//        //  const auto interpro_ids = protein->interproEntries();
//        //  std::string all_ids;
//        //  for (const auto& id : interpro_ids) {
//        //    all_ids = all_ids + id.entry_name + ";";
//        //  }
//        //  report_line = report_line + all_ids + delimiter;
//        //}
//        //else if (token == "%ancestor_domains") {
//        //  const auto ancestors = protein->ancestors();
//        //  std::set<std::string> interpro_terms;
//        //  std::string final_terms;
//        //  for (const auto& ancestor : ancestors) {
//        //    for (const auto& interpro_entry : ancestor.protein->interproEntries()) {
//        //      interpro_terms.insert(interpro_entry.entry_name);
//        //    }
//        //  }  
//
//        //  for (const auto& unique_term : interpro_terms) {
//        //    final_terms = final_terms + unique_term + ";";
//        //  } 
//        //  report_line = report_line + final_terms + delimiter;
//        //}
//        //else if (token == "%ancestor_go") {
//        //  const auto ancestors = protein->ancestors();
//        //  std::set<std::string> go_terms;
//        //  std::string final_terms;
//        //  for (const auto& ancestor : ancestors) {
//        //    for (const auto& go_entry : ancestor.protein->goEntries()) {
//        //      go_terms.insert(go_entry.process_name);
//        //    }
//        //  }  
//
//        //  for (const auto& unique_term : go_terms) {
//        //    final_terms = final_terms + unique_term + ";";
//        //  } 
//        //  report_line = report_line + final_terms + delimiter;
//        //}
//
//        //else if (token == "%ancestor_go_raw") {
//        //  const auto ancestors = protein->ancestors();
//        //  std::set<std::string> go_terms;
//        //  std::string final_terms;
//        //  for (const auto& ancestor : ancestors) {
//        //    for (const auto& go_entry : ancestor.protein->goEntries()) {
//        //      go_terms.insert(go_entry.identifier);
//        //    }
//        //  }  
//
//      //    for (const auto& unique_term : go_terms) {
//      //      final_terms = final_terms + unique_term + ";";
//      //    } 
//      //    report_line = report_line + final_terms + delimiter;
//      //  }
//
//      //  else if (token == "%combined_score") {
//      //    const double cmip_score = protein->score();
//      //    const double instability_index = mipfinder::instability_index(protein->sequence());
//      //    const double combined_score = (1 / (1 + std::exp(-0.1*(instability_index - 70)))) * cmip_score;
//      //    report_line = report_line + std::to_string(combined_score) + delimiter;
//      //  }
//      }
//      of << report_line << "\n";
//    }
//    of.close();
//  }
//}