#ifndef MIPFINDER_HELPERS_H
#define MIPFINDER_HELPERS_H

#include <string>
#include <vector>
#include <algorithm>

namespace mipfinder {
  /* Splits a string into tokens based on the delimiter */
  std::vector<std::string> tokenise(std::string str, char delimiter);

  /* Turns the original string into all-caps */
  void toupper(std::string& str);

  /* Splits a given string at the delimiter into tokens */
  std::vector<std::string> split(const std::string& str, const char delimiter);
}

namespace detail
{
    //TODO: Only really works for containers that have push_back operation but it is fine for
    //this application. Best to wait for std::ranges::to in C++23.
    template <template <typename> typename Container, typename Range>
    requires std::ranges::range<Range>
        Container<std::ranges::range_value_t<Range>> toContainer(Range&& range)
    {
        Container<std::ranges::range_value_t<Range>> container;
        std::ranges::copy(range, std::back_inserter(container));
        return container;
    }
  
}

#endif

