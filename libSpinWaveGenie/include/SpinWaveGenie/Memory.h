#ifndef __SpinWaveGenie_Memory__
#define __SpinWaveGenie_Memory__

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

// implementaiton of std::make_unique for platforms that don't currently have it.
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3656.htm

namespace SpinWaveGenie
{

namespace memory
{

#if __cplusplus == 201402L // C++14
using std::make_unique;

#else // C++11

using namespace std;

template <class T> struct _Unique_if
{
  typedef unique_ptr<T> _Single_object;
};

template <class T> struct _Unique_if<T[]>
{
  typedef unique_ptr<T[]> _Unknown_bound;
};

template <class T, size_t N> struct _Unique_if<T[N]>
{
  typedef void _Known_bound;
};

template <class T, class... Args> typename _Unique_if<T>::_Single_object make_unique(Args &&... args)
{
  return unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <class T> typename _Unique_if<T>::_Unknown_bound make_unique(size_t n)
{
  typedef typename remove_extent<T>::type U;
  return unique_ptr<T>(new U[n]());
}

template <class T, class... Args> typename _Unique_if<T>::_Known_bound make_unique(Args &&...) = delete;

#endif // __cplusplus == 201402L
}
}

#endif /* defined(__SpinWaveGenie_Memory__) */
