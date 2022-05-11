#pragma once

#include <type_traits>
#include <vector>

#include "quantile.hpp"
#include "binning.hpp"



namespace quern
{
	/*
		A histogram over unsigned integer values, from 0 up to some maximum.
			The range is defined by the size of the sample count array.
	*/
	template<
		class T_SampleCounts = std::vector<size_t>>
	struct histogram_by_index
	{
		using sample_counts_t = T_SampleCounts;
		using sample_count_t  = std::decay_t<decltype(std::declval<T_SampleCounts>()[0])>;

		using sample_value_t = bindex_t;
		using quantile_range_t = quantile_range<bindex_t>;


		/*
			Create the histogram (arguments are forwarded to the underlying array)
		*/
		template<typename... Args>
		histogram_by_index(Args&&... args)    : _counts(std::forward<Args>(args)...) {recalculate();}


		size_t         size()       const noexcept    {return _counts.size();}
		sample_count_t population() const noexcept    {return _population;}

		template<typename Index>
		sample_count_t operator[](Index i) const noexcept    {return _counts[i];}

		auto begin()       noexcept    {return _counts.begin();}
		auto begin() const noexcept    {return _counts.begin();}
		auto end()         noexcept    {return _counts.end();}
		auto end()   const noexcept    {return _counts.end();}

		void clear      () noexcept    {_population = 0; for (auto       &c : _counts) c = 0;}
		void recalculate() noexcept    {_population = 0; for (const auto &c : _counts) _population += c;}
		
		void insert(bindex_t insert_index) noexcept
		{
			if (insert_index < 0 || insert_index >= _counts.size()) return;
			++_counts[insert_index];
			++_population;
		}

		void remove(bindex_t remove_index) noexcept
		{
			if (remove_index < 0 || remove_index >= _counts.size()) return;
			--_counts[remove_index];
			--_population;
		}

		void replace(bindex_t insert_index, bindex_t remove_index) noexcept
		{
			insert(insert_index);
			remove(remove_index);
		}


		// Calculate a quantile by scanning the histogram from lowest value to highest.
		quantile_range_t find_quantile(size_t numerator, size_t denominator) const noexcept
		{
			size_t quota = _population * numerator, leq = _counts[0]*denominator;
			size_t index = 0;

			while (index+1 < size() && leq < quota) leq += _counts[++index]*denominator;

			quantile_range_t result;
			result.lower = index;
			if (leq == quota)
				while (index+1 < size()) {if (_counts[++index]) break;}
			result.upper = index;
			return result;
		}

		quantile_range_t find_median() const noexcept    {return find_quantile(1,2);}


	protected:
		sample_counts_t _counts;
		sample_count_t  _population;
	};



	/*
		A histogram over generalized values.
	*/
	template<
		typename T_Sample,
		class T_SampleCounts = std::vector<size_t>>
	class histogram
	{
	public:
		using sample_counts_t = T_SampleCounts;
		using sample_count_t  = std::decay_t<decltype(std::declval<T_SampleCounts>()[0])>;

		using data_t = histogram_by_index<T_SampleCounts>;

		using sample_value_t = T_Sample;
		using quantile_range_t = quantile_range<sample_value_t>;


	public:
		/*
			Create the histogram (arguments are forwarded to the underlying array)
		*/
		template<typename... Args>
		histogram(Args&&... args)    : data(std::forward<Args>(args)...) {recalculate();}


		size_t         size()       const noexcept    {return data.size();}
		sample_count_t population() const noexcept    {return data.population();}


		// TODO map-like iteration


		void clear      () noexcept    {data.clear();}
		void recalculate() noexcept    {data.recalculate();}

		
		void insert(sample_value_t new_v) noexcept
		{
			
			++_counts[insert_index];
			++_population;
		}

		void remove(sample_value_t old_v) noexcept
		{
			--_counts[remove_index];
			--_population;
		}

		void replace(sample_value_t new_v, sample_value_t old_v) noexcept
		{
			insert(insert_index);
			remove(remove_index);
		}


		// Calculate a quantile by scanning the histogram from lowest value to highest.
		quantile_range_t find_quantile(size_t numerator, size_t denominator) const noexcept
		{
			size_t quota = _population * numerator, leq = _counts[0]*denominator;
			size_t index = 0;

			while (index+1 < size() && leq < quota) leq += _counts[++index]*denominator;

			quantile_range_t result;
			result.lower = index;
			if (leq == quota)
				while (index+1 < size()) {if (_counts[++index]) break;}
			result.upper = index;
			return result;
		}

		quantile_range_t find_median() const noexcept    {return find_quantile(1,2);}


	public:
		data_t data;
	};
}
