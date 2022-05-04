#pragma once

#include <exception>
#include <utility>
#include <type_traits>
#include <vector>



namespace edb
{
	namespace detail
	{
		template<typename T> using std_vector_of = std::vector<T>;
	}

	/*
		Represents the location of a quantile.
			When samples are evenly divided, this can be an exclusive range
			containing no samples (such as the space between two histogram slots).
	*/
	template<typename T_Value>
	struct quantile_range
	{
		using value_t = T_Value;

		value_t lower, upper;

		bool is_range() const noexcept    {return lower != upper;}
		bool is_value() const noexcept    {return lower == upper;}

		operator      float () const noexcept    {return .5f*(lower+upper);}
		operator      double() const noexcept    {return .5 *(lower+upper);}
		operator long double() const noexcept    {return .5L*(lower+upper);}
	};

	/*
		Represents a fraction defining a quantile.
	*/
	template<typename T_Int = size_t>
	struct quantile_fraction
	{
		using int_t = T_Int;

		quantile_fraction(int_t numerator, int_t denominator = 1)
			:
			num(numerator), den(denominator) {}

		union
		{
			struct {int_t num,       den;};
			struct {int_t numerator, denominator;};
		};

		bool operator< (const quantile_fraction &o) const noexcept    {return num*o.den <  o.num*den;}
		bool operator<=(const quantile_fraction &o) const noexcept    {return num*o.den <= o.num*den;}
		bool operator> (const quantile_fraction &o) const noexcept    {return num*o.den >  o.num*den;}
		bool operator>=(const quantile_fraction &o) const noexcept    {return num*o.den >= o.num*den;}
		bool operator==(const quantile_fraction &o) const noexcept    {return num*o.den == o.num*den;}
		bool operator!=(const quantile_fraction &o) const noexcept    {return num*o.den != o.num*den;}

		operator      float () const noexcept    {return num / float(den);}
		operator      double() const noexcept    {return num / double(den);}
		operator long double() const noexcept    {return num / long double(den);}
	};


	/*
		A simple histogram for non-negative integer values.
			Essentially an array of non-negative sample counts.
	*/
	template<class T_SampleCounts = std::vector<size_t>>
	struct histogram_basic
	{
		using sample_counts_t = T_SampleCounts;
		using sample_count_t  = std::decay_t<decltype(std::declval<T_SampleCounts>()[0])>;


		using quantile_range_t = quantile_range<size_t>;


		histogram_basic()    : _population(0) {for (auto &c : _counts) c = 0;}


		size_t         size()  const noexcept    {return _counts.size();}
		sample_count_t population() const noexcept    {return _population;}


		template<typename Index>
		sample_count_t operator[](Index i) const noexcept    {return _counts[i];}

		
		auto begin()       noexcept    {return _counts.begin();}
		auto begin() const noexcept    {return _counts.begin();}
		auto end()         noexcept    {return _counts.end();}
		auto end()   const noexcept    {return _counts.end();}

		
		template<typename Index>
		void insert(Index insert_index) noexcept
		{
			if (insert_index < 0 || insert_index >= _counts.size()) return;
			++_counts[insert_index];
			++_population;
		}

		template<typename Index>
		void remove(Index remove_index) noexcept
		{
			if (remove_index < 0 || remove_index >= _counts.size()) return;
			--_counts[remove_index];
			--_population;
		}

		template<typename Index>
		void replace(Index insert_index, Index remove_index) noexcept
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
		An object that tracks a quantile in a histogram.
	*/
	template<typename T_Index = size_t, typename T_SampleCount = size_t>
	struct quantile_position
	{
		using index_t = T_Index;
		using sample_count_t = T_SampleCount;


		// Definition of the quantile.
		//  eg, with denominator 100 (percentile), numerator 50 tracks the median.
		quantile_fraction<T_Index> quantile;

		// The lower and upper bins of the quantile.  lower <= upper.
		//   These may differ if the histogram has an even number of samples
		//   and those samples are evenly divided between two sub-ranges.
		quantile_range<index_t> range;


		// This value tracks how many samples are below bin_upper.
		sample_count_t samples_lower;


		// DEBUG
		short last_adjust = 0;


	public:
		/*
			Recalculate this quantile based on the given histogram.
				Re-calculates samples_lower and calls adjust.
		*/
		template<typename Histogram>
		void recalculate(const Histogram &h, index_t hint_index = 0)
		{
			if (hint_index >= h.size()) hint_index = h.size() - (h.size() > 0);

			range.lower = range.upper = hint_index;
			samples_lower = 0;
			for (index_t i = 0; i < hint_index; ++i) samples_lower += h[i];
			adjust(h);
		}

		/*
			Adjust this quantile based on the given histogram.
				Unlike recalculate() this assumes samples_lower is kept up-to-date.
		*/
		template<typename Histogram>
		void adjust(const Histogram &h)
		{
			// "smash" any range to its upper bound
			index_t bin = range.upper;
			sample_count_t population = h.population();
			sample_count_t
				here  = h[bin];
			size_t
				gte   = population - samples_lower,
				lte   = here + samples_lower;
			size_t
				lte_ratio = population*quantile.num,
				gte_ratio = population*(quantile.den-quantile.num);

			if (lte*quantile.den < lte_ratio)
			{
				last_adjust = 1;

				// Slide the quantile higher
				while (bin+1 < h.size() && lte*quantile.den < lte_ratio)
				{
					samples_lower += here;
					here = h[++bin];
					//q.samples_higher -= here;
					lte += here;
					if (lte*quantile.den >= lte_ratio) break;
				}

				// Determine the median bin, or bin range in case of a split
				range.lower = bin;
				if (lte*quantile.den == lte_ratio)
				{
					samples_lower += here;
					while (bin+1 < h.size() && h[++bin] == 0) {}
				}
				range.upper = bin;
			}
			else if (gte*quantile.den < gte_ratio)
			{
				last_adjust = -1;

				// Slide the quantile lower
				while (bin > 0 && gte*quantile.den < gte_ratio)
				{
					//q.samples_higher += here;
					here = h[--bin];
					samples_lower -= here;
					gte += here;
					if (gte*quantile.den >= gte_ratio) break;
				}

				// Determine the median bin, or bin range in case of a split
				range.upper = bin;
				if (gte*quantile.den == gte_ratio)
				{
					while (bin > 0 && h[--bin] == 0) {}
				}
				range.lower = bin;
			}
			else
			{
				last_adjust = 0;

				// Elements <= bin and >= bin are partitioned...
				range.lower = range.upper = bin;
				
				while (range.lower > 0) // expand range downward
				{
					lte -= h[range.lower];
					if (lte*quantile.den < lte_ratio) break;
					--range.lower;
				}
				while (range.upper+1 < h.size()) // expand range upward
				{
					gte -= h[range.upper];
					if (gte*quantile.den < gte_ratio) break;
					samples_lower += h[range.upper];
					++range.upper;
				}
			}
		}
	};


	/*
		A histogram over non-negative integers which tracks the values
			of various quantiles with each update.
		
		array_t          -- used to store histogram counts
		sample_count_t   -- used to store histogram counts
		index_t          -- must be able to store twice the number of bins in the histogram
		quantile_array_t -- type used to store quantile data
	*/
	template<
		class    T_HistogramBase = histogram_basic<>,
		class    T_Quantiles     = std::vector<quantile_position<size_t, typename T_HistogramBase::sample_count_t>>>
	class quantile_tracker
	{
	public:
		using histogram_t        = T_HistogramBase;
		using sample_count_t     = typename T_HistogramBase::sample_count_t;
		using quantiles_t        = T_Quantiles;
		using quantile_t         = std::decay_t<decltype(std::declval<T_Quantiles>()[0])>;
		using index_t            = typename quantile_t::index_t;


	public:
		/*
			Create the histogram and calculate initial quantiles.
				The list of quantiles must be sorted.
		*/
		quantile_tracker(quantiles_t quantiles)
			:
			_quantiles(std::move(quantiles))
		{
			// TODO calculate initial quantiles.
			for (auto &q : _quantiles)
			{
				if (q.quantile.den == 0)             throw std::logic_error("Invalid quantile: denominator = 0");
				if (q.quantile.num == 0)             throw std::logic_error("Invalid quantile: numerator = 0");
				if (q.quantile.num > q.quantile.den) throw std::logic_error("Invalid quantile: numerator > denominator");

				q.recalculate(_histogram);
			}
		}


		/*
			Access histogram and quantile readouts.
		*/
		const histogram_t &histogram() const noexcept    {return _histogram;}
		const quantiles_t &quantiles() const noexcept    {return _quantiles;}


		/*
			Insert an item.
		*/
		void insert(index_t insert_index)
		{
			_histogram.insert(insert_index);
			for (auto &q : _quantiles)
			{
				if (insert_index < q.range.upper) ++q.samples_lower;
				q.adjust(_histogram);
			}
		}
		void remove(index_t remove_index)
		{
			_histogram.remove(remove_index);
			for (auto &q : _quantiles)
			{
				if (remove_index < q.range.upper) --q.samples_lower;
				q.adjust(_histogram);
			}
		}

		/*
			Replace an item.
				Essentially "moves" a sample to the insert index from the remove index.
				This can save work for quantiles that don't need updating.
		*/
		void replace(index_t insert_index, index_t remove_index)
		{
			if (insert_index != remove_index)
			{
				_histogram.replace(insert_index, remove_index);

				for (auto &q : _quantiles)
				{
					//auto // we can skip if the samples are on the same side of the quantile
					//	i_chg = (insert_index > q.range.lower) - (insert_index < q.range.upper),
					//	r_chg = (remove_index > q.range.lower) - (remove_index < q.range.upper);

					//if (i_chg == r_chg) continue; // TODO fix this

					// Adjust the quantile.
					q.samples_lower += (insert_index < q.range.upper) - (remove_index < q.range.upper);
					q.adjust(_histogram);
				}
			}
		}


	private:
		histogram_t    _histogram;
		index_t        _q_denominator;
		quantiles_t    _quantiles;
	};
}
