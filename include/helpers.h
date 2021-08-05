#ifndef MIPFINDER_HELPERS_H
#define MIPFINDER_HELPERS_H

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

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
    template <typename T>
    concept PushBackableContainer = requires (T t, typename T::value_type v)
    {
        { t.push_back(v) } -> std::same_as<void>;
    };

    template <typename T>
    concept InsertableContainer = requires (T t, typename T::value_type v)
    {
        { t.insert(v) } -> std::same_as<std::pair<typename T::iterator, bool>>;
    };

    template <typename T>
    concept EmplaceableContainer = requires (T t, typename T::value_type v)
    {
        t.emplace_back(v);
    };

    //TODO: Only really works for containers that have push_back operation but it is fine for
    //this application. Best to wait for std::ranges::to in C++23.
    template <template <typename> typename Container, typename Range>
    requires std::ranges::range<Range>
        Container<std::ranges::range_value_t<Range>> toContainer(Range&& range)
    {
        Container<std::ranges::range_value_t<Range>> container;
        if constexpr (PushBackableContainer<Range>) {
            std::ranges::copy(range, std::back_inserter(container));
        }
        else if constexpr (InsertableContainer<Range>) {
            std::ranges::copy(range, std::inserter(container));
        }
        return container;
    }
  
}

#endif

