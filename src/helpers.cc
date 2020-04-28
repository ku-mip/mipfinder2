#include <sstream>

#include "helpers.h"

namespace mipfinder {
  std::vector<std::string> tokenise(std::string str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream s{str};
    std::string string_token{};
    while(std::getline(s, string_token, delimiter)) {
      if (string_token.length() != 0) {
        tokens.emplace_back(string_token);
      }
    }
    return tokens;
  }

  void toupper(std::string& str) {
    for (char& c :str) {
      c = std::toupper(c);
    }
  }

  std::vector<std::string> split(const std::string& str, const char delimiter)
  {
    std::vector<std::string> split_strings;
    for (std::size_t current_pos = 0;;) {
      const auto token_pos = str.find(delimiter, current_pos);
      if (token_pos == std::string::npos) {
        break;
      }

      const auto next_token_pos = str.find(delimiter, token_pos + 1);
      const auto substr_len = next_token_pos - token_pos;
      split_strings.push_back(str.substr(token_pos, substr_len));
      current_pos = next_token_pos;
    }
    return split_strings;
  }
}

