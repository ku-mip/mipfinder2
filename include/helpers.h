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
    template <typename T, typename U>
    requires std::ranges::range<T>
    auto toContainer(T&& range)
    {
        using View = decltype(std::declval<T>().base());
        using ContainerType = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<View>().base())>>;
        ContainerType container;
        // auto container_size = std::distance(std::begin(range_view), std::end(range_view));
        // container.resize(container_size);
        std::ranges::copy(range_view, std::back_inserter(container));
        return container;
    }
  
}

#endif

