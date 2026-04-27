/* Copyright (c) 2015 Tobias Becker
**
**
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is
** furnished to do so, subject to the following conditions:
**
**
**
** The above copyright notice and this permission notice shall be included in
** all copies or substantial portions of the Software.
**
** Adapted by Quentin Schau on 23/04/2026 to work with pair of double
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
** THE SOFTWARE.
*/

#pragma once

#include <vector>
#include <array>
#include <limits>
#include <type_traits>
#include <utility>
#include <numeric>
#include <algorithm>
#include <cstddef>

/* PROPERTIES:
** - Inplace
** - O(k * n)
** - O(n) temporary memory needed.
** - Stable sort due to the fact that it is a LSB radix sort.
** - Works with forward iterators but takes advantage of random access iterators
**   with contiguous memory to save unnecessary copies of the temporary array.
** - Contiguous memory mode is deduced from the iterator type being random access.
**   If your data struct does not fit in, specify manually! Otherwise, at least in
**   debug mode an assertion should be triggered.
** - All byte-sized histograms are generated in scan over the input.
** - A fully sorted input is detected in this first pass and sorting is done.
** - During histogram generation, pathological rounds are remembered to be skipped
**   later because they wouldn't change any ordering anyway. For instance sorting n
**   unsigned integers that are all < 256 would only take O(2 * n) instead of O(5 * n).
** - All three cases (unsigned integer, signed integer and (signed) floating point)
**   are handled in one block without runtime overhead that stems from unused
**   cases. The compiler should optimize unneeded checks away because they're known
**   at compile time.
** - Floating point bit manipulation is compliant with C++'s strict alias rules.
** - Can switch to insertion sort for small n.
*/

//#define INSERTION_SORT_FOR_SMALL_N

namespace RadixSortDetails
{
	using byte = unsigned char;

	template<class T>
	using lim = typename std::numeric_limits<T>;

	template<class Iter>
	using iterator_value_type = typename std::iterator_traits<Iter>::value_type;

	enum
	{
		HISTOGRAM_SIZE = lim<byte>::max() + 1,
		HISTOGRAM_HALF = HISTOGRAM_SIZE / 2
	};

	template<class T>
	/// Extract byte from any integer (shift = 0 for LSB).
	inline byte const * ubyte(T const & num)
	{
		static_assert(std::is_arithmetic<T>::value, "Only for fundamental number types that are allowed.");
		static_assert(lim<T>::is_integer || lim<T>::is_iec559,
			"Only know how to handle ints and (32-bit and 64-bit) IEEE-754 floating point values.");

		/*
		 * This is in accordance with the standard, which allowes casting any pointer to an unsigned
		 * char pointer for inspection (i.e. reading). Using memcpy would have been the other option,
		 * but as we only need one byte anyway and this is easily available though indexing this seemed
		 * like the better choice.
		 */
		return reinterpret_cast<byte const *>(&num);
	}

	template<class Iter>
	bool generate_histograms(Iter const begin, Iter const end, typename std::iterator_traits<Iter>::difference_type const n,
		std::array<std::array<int, HISTOGRAM_SIZE>, sizeof(iterator_value_type<Iter>)> & histograms)
	{
		assert(histograms.data());

		bool sorted = true;
		int max_counts[sizeof(iterator_value_type<Iter>)]{};

		auto previous = begin;
		// One-time pass over input to build all histograms at once. Also check if sorting is necessary.
		for (auto iter = begin; iter != end; ++iter)
		{
			if (*previous > *iter)
			{
				sorted = false;
			}
			previous = iter;

			auto const positions = ubyte(*iter);
			int shift = 0;
			for (auto & hist : histograms)
			{
				if (++hist[positions[shift]] > max_counts[shift])
				{
					max_counts[shift] = hist[positions[shift]];
				}
				++shift;
			}
		}

		for (auto const & h : histograms)
		{
			assert(std::accumulate(h.begin(), h.end(), 0) == n);
		}

		for (int i = 0; i < sizeof(iterator_value_type<Iter>); ++i)
		{
			if (max_counts[i] == n)
			{
				histograms[i][0] = -1;
			}
		}

		return sorted;
	}

	template<class Iter, class Iter2>
	/// For non-contiguous memory (last parameter), we always have to copy in between rounds.
	inline void post_round_copy(Iter & begin, Iter,
		Iter2 & begin2, Iter2 & end2,
		bool, std::false_type)
	{
		std::copy(begin2, end2, begin);
	}

	template<class Iter>
	/// For contiguous memory (last parameter), we can replace in-between copying by swapping iterators.
	inline void post_round_copy(Iter & begin, Iter & end,
		Iter & begin2, Iter & end2,
		bool & swapped, std::true_type)
	{
		using namespace std;

		/*
		 * We can only do this check here, because (forward_)list will never get here and we *know* that one
		 * of the ranges comes from a vector, which is contiguous.
		 *
		 * So this can only fail for some random access data structure that does not use contiguous memory.
		 * But then again the template parameter CONTIGUOUS should have been manually set to false and we
		 * would not have come here.
		 */
		assert(distance(begin, end) > 1
			&& distance(begin, end) == &end[-1] - &*begin + 1
			&& &end[-1] - &*begin == &end2[-1] - &*begin2);

		swap(begin, begin2);
		swap(end, end2);
		swapped = !swapped;
	}

#ifdef INSERTION_SORT_FOR_SMALL_N
	template<class Iter>
	void insertionsort(Iter const begin, Iter const end)
	{
		using namespace std;

		if (begin == end)
		{
			return;
		}

		for (auto unsorted = next(begin); unsorted != end; ++unsorted)
		{
			auto const current = *unsorted;
			auto shift = unsorted;

			do
			{
				auto before = prev(shift);
				if (*before > current)
				{
					*shift = *before;
					shift = move(before);
				}
				else
				{
					break;
				}
			} while (shift != begin);

			*shift = current;
		}
	}
#endif
}

/// Radix sort for all (un-)signed integer types and 32-bit float and 64-bit double (as long as IEEE-754 compatible).
template<class Iter,
	bool CONTIGUOUS = std::is_same<typename std::iterator_traits<Iter>::iterator_category,
	std::random_access_iterator_tag>::value>
	void radixsort(Iter begin, Iter end)
{
	using namespace std;
	using namespace RadixSortDetails;

	using T = typename iterator_traits<Iter>::value_type;
	using limT = lim<T>;

	static_assert(sizeof(T) <= 8, "Not made for an extended long double.");

	auto const n = distance(begin, end);

#ifdef INSERTION_SORT_FOR_SMALL_N
	if (n < 48 * sizeof(T))
	{
		insertionsort(begin, end);
		return;
	}
#else
	if (n < 2)
	{
		return;
	}
#endif

	// Holds a handfull of histograms.
	array<array<int, HISTOGRAM_SIZE>, sizeof(T)> histograms{};

	// This returns true, if input is already sorted. Then we don't have to do anything.
	if (!generate_histograms(begin, end, n, histograms))
	{
		// Count the number of negative HSBs.
		auto const & last = histograms.back();
		auto const negatives = limT::is_signed
			? accumulate(last.cbegin() + HISTOGRAM_HALF, last.cend(), 0)
			: 0;

		int offsets[HISTOGRAM_SIZE];
		vector<T> tmp(n);
		auto begin2 = tmp.begin();

		// These two are only needed for CONTIGUOUS
		auto end2 = tmp.end();
		bool swapped = false;

		int shift_bytes = 0;
		for (auto const & hist : histograms)
		{
			if (hist.front() >= 0)
			{
				// If we need special care, it will only be in the last round...
				bool const last_round = &hist[0] == &last[0];
				// ... and only if we have a signed type.
				bool const special = last_round && limT::is_signed;

				// Negative values need come first. We do this by shifting the first positive offset.
				*offsets = special ? negatives : 0;

				// The first (positive half is the same for everyone.
				int i;
				for (i = 1; i < HISTOGRAM_HALF; ++i)
				{
					offsets[i] = offsets[i - 1] + hist[i - 1];
				}

				if (special && limT::is_iec559)
				{
					/* This is for floating point values. Without special care, they'd come out
					 * sorted by absolute value after the positive ones, like this:
					 *
					 * 1, 3, 5, -2, -4, -6
					 *
					 * So the 0 offset goes to the end and we start summing backwards (but only for
					 * the second (negative) half).
					 */

					// Downward upper half...
					offsets[HISTOGRAM_SIZE - 1] = 0;
					for (i = HISTOGRAM_SIZE - 2; i >= HISTOGRAM_HALF; --i)
					{
						offsets[i] = offsets[i + 1] + hist[i + 1];
					}
					// ... and upward lower half.
					for (i = HISTOGRAM_HALF; i < HISTOGRAM_SIZE; ++i)
					{
						offsets[i] += hist[i];
					}
				}
				else
				{
					// For signed integers, the offsets can be build in increasing order as for unsigned ones.
					if (special)
					{
						// Except the signed integers need the 09 offset at the beginning of the second (negative) half.
						offsets[i++] = 0;
					}

					// After that they're the same again.
					for (; i < HISTOGRAM_SIZE; ++i)
					{
						offsets[i] = offsets[i - 1] + hist[i - 1];
					}
				}

				// Now values are written to the offsets previously computed.
				// This establishes the order, byte by byte.
				for (auto iter = begin; iter != end; ++iter)
				{
					auto const radix = ubyte(*iter)[shift_bytes];
					if (special && limT::is_iec559 && radix >= HISTOGRAM_HALF)
					{
						begin2[--offsets[radix]] = *iter;
					}
					else
					{
						begin2[offsets[radix]++] = *iter;
					}
				}

				// Either swap iterator roles for CONTIGUOUS, or copy from tmp to input...
				post_round_copy(begin, end, begin2, end2, swapped, integral_constant<bool, CONTIGUOUS>());
			}

			// ... and repeat with the next higher byte.
			++shift_bytes;
		}

		// We might end up with swapped roles. Then we have to do an additional copy.
		if (CONTIGUOUS && swapped)
		{
			// Temporary vector will die after this anyway.
			move(begin, end, begin2);
		}
	}
}

template<class C>
/// Sort given container range with radix sort.
inline C & radixsort(C && cont)
{
	using namespace std;
	radixsort(begin(cont), end(cont));
	return cont;
}

namespace RadixSortByFirstDetails
{
    using byte = unsigned char;

    template<class T>
    using lim = std::numeric_limits<T>;

    template<class T>
    inline const byte* ubyte(const T& num)
    {
        static_assert(std::is_arithmetic_v<T>,
            "radixsort_by_first requires an arithmetic key type.");
        static_assert(lim<T>::is_integer || lim<T>::is_iec559,
            "Only integer and IEEE-754 float/double keys are supported.");

        return reinterpret_cast<const byte*>(&num);
    }
}

template<class Pair>
void radixsort_by_first(std::vector<Pair>& v)
{
    using namespace RadixSortByFirstDetails;

    using Key = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<Pair>().first)>>;

    static_assert(std::is_arithmetic_v<Key>,
        "Pair::first must be an arithmetic type.");
    static_assert(sizeof(Key) <= 8,
        "This implementation is intended for up to 64-bit keys.");

    constexpr std::size_t HISTOGRAM_SIZE = 256;
    constexpr std::size_t HISTOGRAM_HALF = 128;
    constexpr std::size_t KEY_BYTES = sizeof(Key);

    const std::size_t n = v.size();
    if (n < 2) return;

    std::array<std::array<std::size_t, HISTOGRAM_SIZE>, KEY_BYTES> histograms{};
    std::array<std::size_t, KEY_BYTES> max_counts{};

    bool sorted = true;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (i > 0 && v[i - 1].first > v[i].first)
            sorted = false;

        const byte* positions = ubyte(v[i].first);
        for (std::size_t b = 0; b < KEY_BYTES; ++b)
        {
            const std::size_t c = ++histograms[b][positions[b]];
            if (c > max_counts[b]) max_counts[b] = c;
        }
    }

    if (sorted) return;

    constexpr std::size_t SKIP_ROUND = std::numeric_limits<std::size_t>::max();
    for (std::size_t b = 0; b < KEY_BYTES; ++b)
    {
        if (max_counts[b] == n)
            histograms[b][0] = SKIP_ROUND;
    }

    const auto& last = histograms[KEY_BYTES - 1];
    const std::size_t negatives =
        lim<Key>::is_signed
            ? std::accumulate(last.begin() + HISTOGRAM_HALF, last.end(), std::size_t{0})
            : 0;

    std::vector<Pair> tmp(n);
    auto* src = &v;
    auto* dst = &tmp;

    std::array<std::size_t, HISTOGRAM_SIZE> offsets{};

    for (std::size_t shift_bytes = 0; shift_bytes < KEY_BYTES; ++shift_bytes)
    {
        const auto& hist = histograms[shift_bytes];
        if (hist[0] == SKIP_ROUND) continue;

        const bool last_round = (shift_bytes == KEY_BYTES - 1);
        const bool special = last_round && lim<Key>::is_signed;

        offsets[0] = special ? negatives : 0;

        std::size_t i = 1;
        for (; i < HISTOGRAM_HALF; ++i)
            offsets[i] = offsets[i - 1] + hist[i - 1];

        if (special && lim<Key>::is_iec559)
        {
            offsets[HISTOGRAM_SIZE - 1] = 0;
            for (int j = static_cast<int>(HISTOGRAM_SIZE) - 2; j >= static_cast<int>(HISTOGRAM_HALF); --j)
                offsets[j] = offsets[j + 1] + hist[j + 1];

            for (std::size_t j = HISTOGRAM_HALF; j < HISTOGRAM_SIZE; ++j)
                offsets[j] += hist[j];
        }
        else
        {
            if (special)
                offsets[i++] = 0;

            for (; i < HISTOGRAM_SIZE; ++i)
                offsets[i] = offsets[i - 1] + hist[i - 1];
        }

        for (const auto& elem : *src)
        {
            const byte radix = ubyte(elem.first)[shift_bytes];

            if (special && lim<Key>::is_iec559 && radix >= HISTOGRAM_HALF)
                (*dst)[--offsets[radix]] = elem;
            else
                (*dst)[offsets[radix]++] = elem;
        }

        std::swap(src, dst);
    }

    if (src != &v)
        v = std::move(*src);
}