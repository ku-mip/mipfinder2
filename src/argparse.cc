#include <vector>
#include "argparse.h"

namespace mipfinder
{
    std::vector<std::string> split(std::string str, std::string delimiter)
    {
        std::vector<std::string> tokens;
        std::string::size_type pos = 0;
        std::string::size_type prev_pos = 0;
        while ((pos = str.find(delimiter, pos)) != std::string::npos) {
            auto token = str.substr(prev_pos, pos - prev_pos);
            pos += delimiter.length();
            prev_pos = pos;

            if (token.empty()) {
                continue;
            }
            tokens.push_back(token);
        }
        auto final_token = str.substr(prev_pos, std::string::npos);
        tokens.push_back(final_token);
        return tokens;
    }

    Args parse_args(int argc, char* argv[])
    {
        std::string arguments;
        for (int i = 0; i < argc; i++) {
            arguments.append(argv[i]).append(" ");
        }
        return parse_args(arguments);
    }



    /* A very simple command line parser that extracts options. Every option is preceded by a hyphen ('-')
     * token followed by an optional value, e.g. "-i file.txt".
     */
    Args parse_args(const std::string& args)
    {
        Args parsed_args;
        auto tokens = split(args, " ");
        for (auto it = std::begin(tokens); it != std::end(tokens); ) {
            auto& token = it;
            auto next_token = std::next(it);

            if (next_token == std::end(tokens)) {
                parsed_args.insert({*token, ""});
                break;
            }

            if ((*token).front() == '-') { //handles both short and long versions (e.g. -h vs --help)
                if ((*next_token).front() != '-') {
                    parsed_args.insert({*token, *next_token});
                    std::advance(it, 2);
                    continue;
                }
            }
            parsed_args.insert({*token, ""});
            std::advance(it, 1);
        }
        return parsed_args;
    }
}