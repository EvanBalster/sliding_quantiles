#pragma once

#include <exception>
#include <utility>
#include <array>

#include "quantile.hpp"
#include "histogram.hpp"



namespace quern
{
	namespace detail
	{
		template<typename T> using std_vector_of = std::vector<T>;
	}


	/*
		A histogram which tracks total population and quantile values.
		
		array_t          -- used to store histogram counts
		sample_count_t   -- used to store histogram counts
		index_t          -- must be able to store twice the number of bins in the histogram
		quantile_array_t -- type used to store quantile data
	*/
	template<
		class    T_HistogramBase>
		//class    T_Quantiles      = std::vector<quantile_fraction<typename T_HistogramBase::index_t>>,
		//typename T_QuantileValues = std::vector<tracked_quantile<typename T_HistogramBase::count_t, typename T_HistogramBase::index_t>>>
	class histogram_tracked
	{
	public:
		static_assert(T_HistogramBase::dimensionality == 1, "histogram_tracked must be 1-dimensional");

		using histogram_t  = T_HistogramBase;
		using sample_t     = typename histogram_t::sample_t;
		using count_t      = typename histogram_t::count_t;
		using index_t      = typename histogram_t::index_t;
		using binning_t    = typename histogram_t::binning_t;
		using params_t     = typename histogram_t::params_t;

		/*
			Data structure representing one tracked quantile.
		*/
		struct quantile
		{
			// Definition of the quantile.
			quantile_fraction<index_t> quantile;

			// The lower and upper bins of the quantile.  lower <= upper.
			//   These may differ if the histogram has an even number of samples
			//   and those samples are evenly divided between two sub-ranges.
			quantile_range<index_t> index_range;


			// This value tracks how many samples are below bin_upper.
			count_t samples_lower;


			// DEBUG
			short last_adjust = 0;


			void recalculate(const histogram_t &h, count_t population, bindex_t hint_index = 0);
			void adjust     (const histogram_t &h, count_t population);
		};

		using quantiles_t = std::vector<quantile>;


	public:
		/*
			Default constructor.  This empty histogram will not accept samples.
		*/
		explicit histogram_tracked()    : _histogram() {}

		/*
			Set up empty bins based on an array of binning rules.
		*/
		histogram_tracked(const binning_t &binning)    : _histogram(binning), _population(0) {}
		histogram_tracked(const params_t  &params )    : _histogram(params ), _population(0) {}

		/*
			As above but also specify quantiles to track in the constructor.
		*/
		template<typename QuantileList>
		histogram_tracked(const binning_t &binning, const QuantileList &quantiles)    : _histogram(binning, count_t(0)), _population(0) {_init_quantiles(quantiles);}
		template<typename QuantileList>
		histogram_tracked(const params_t  &params , const QuantileList &quantiles)    : _histogram(params , count_t(0)), _population(0) {_init_quantiles(quantiles);}


		template<typename QuantileList>
		void add_quantiles(const QuantileList &quantiles)
		{
			_quantiles.reserve(_quantiles.size()+std::size(quantiles));
			for (auto &q : quantiles)
			{
				_quantiles.emplace_back(quantile{q});
				_quantiles.back().recalculate(_histogram, _population);
			}
		}

		void recalculate()
		{
			_population = _histogram.calc_population();

			_quantiles.resize(_quantiles.size());

			for (auto &q : _quantiles)
			{
				q.recalculate(_histogram, _population);
			}
		}


		/*
			Access histogram and quantile readouts.
		*/
		const histogram_t &histogram() const noexcept    {return _histogram;}
		const quantiles_t &quantiles() const noexcept    {return _quantiles;}


		const count_t     population() const noexcept    {return _population;}


		/*
			Insert an item.
		*/
		void insert(sample_t new_sample)
		{
			count_t miss = 0;
			_histogram.at(new_sample, miss) += 1;
			if (!miss)
			{
				++_population;
				for (auto &q : _quantiles)
				{
					if (new_sample < q.index_range.upper) ++q.samples_lower;
					q.adjust(_histogram, _population);
				}
			}
			else {for (auto &q : _quantiles) q.last_adjust = -2;}
		}
		void remove(sample_t old_sample)
		{
			count_t hit = 1;
			_histogram.at(old_sample, hit) -= 1;
			if (hit)
			{
				--_population;
				for (auto &q : _quantiles)
				{
					if (old_sample < q.index_range.upper) --q.samples_lower;
					q.adjust(_histogram, _population);
				}
			}
			else {for (auto &q : _quantiles) q.last_adjust = -3;}
		}

		/*
			Replace an item.
				Essentially "moves" a sample to the insert index from the remove index.
				This can save work for quantiles that don't need updating.
		*/
		void replace(sample_t new_sample, index_t old_sample)
		{
			index_t new_index = _histogram.index_for(new_sample);
			if (new_index == BIN_REJECT)
			{
				remove(old_sample);
				return;
			}
			index_t old_index = _histogram.index_for(old_sample);
			if (old_index == BIN_REJECT)
			{
				insert(new_sample);
				return;
			}

			if (new_index != old_index)
			{
				count_t dummy;
				_histogram.at_index(new_sample, dummy) += 1;
				_histogram.at_index(old_sample, dummy) -= 1;

				for (auto &q : _quantiles)
				{
					q.last_adjust = 9;

					// No need to adjust if samples are both outside the quantile in the same direction
					if (new_index > q.index_range.upper && old_index > q.index_range.upper) continue;
					if (new_index < q.index_range.lower && old_index < q.index_range.lower) continue;

					// Adjust the quantile.
					q.samples_lower += (new_sample < q.index_range.upper) - (old_sample < q.index_range.upper);
					q.adjust(_histogram, _population);
				}
			}
		}


	private:
		template<typename QuantileList>
		void _init_quantiles(const QuantileList &quantiles)
		{
			_quantiles.reserve(std::size(quantiles));
			for (auto &q : quantiles) _quantiles.emplace_back(quantile{q, {0,_histogram.bins()-1}});
		}

		histogram_t    _histogram;
		count_t        _population;
		quantiles_t    _quantiles;
	};
}



template<typename Histogram>
void quern::histogram_tracked<Histogram>::quantile::recalculate
	(const Histogram &h, count_t population, bindex_t hint_index)
{
	if (quantile.den <= 0)            throw std::logic_error("Invalid quantile: denominator <= 0");
	if (quantile.num <= 0)            throw std::logic_error("Invalid quantile: ratio <= 0");
	if (quantile.num >= quantile.den) throw std::logic_error("Invalid quantile: ratio >= 1");

	auto size = h.bins();

	if (hint_index >= size) hint_index = size - (size > 0);

	index_range.lower = index_range.upper = hint_index;
	samples_lower = 0;
	for (index_t i = 0; i < hint_index; ++i) samples_lower += h.count_at(i);
	adjust(h, population);
}

template<typename Histogram>
void quern::histogram_tracked<Histogram>::quantile::adjust
	(const histogram_t &h, count_t population)
{
	auto size = h.bins();

	// "smash" any range to its upper bound
	bindex_t bin = index_range.upper;
	count_t
		here  = h.count_at(bin);
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
		while (bin+1 < size && lte*quantile.den < lte_ratio)
		{
			samples_lower += here;
			here = h.count_at(++bin);
			//q.samples_higher -= here;
			lte += here;
			if (lte*quantile.den >= lte_ratio) break;
		}

		// Determine the median bin, or bin range in case of a split
		index_range.lower = bin;
		if (lte*quantile.den == lte_ratio)
		{
			samples_lower += here;
			while (bin+1 < size && h.count_at(++bin) == 0) {}
		}
		index_range.upper = bin;
	}
	else if (gte*quantile.den < gte_ratio)
	{
		last_adjust = -1;

		// Slide the quantile lower
		while (bin > 0 && gte*quantile.den < gte_ratio)
		{
			//q.samples_higher += here;
			here = h.count_at(--bin);
			samples_lower -= here;
			gte += here;
			if (gte*quantile.den >= gte_ratio) break;
		}

		// Determine the median bin, or bin range in case of a split
		index_range.upper = bin;
		if (gte*quantile.den == gte_ratio)
		{
			while (bin > 0 && h.count_at(--bin) == 0) {}
		}
		index_range.lower = bin;
	}
	else
	{
		last_adjust = 0;

		// Elements <= bin and >= bin are partitioned...
		index_range.lower = index_range.upper = bin;
				
		while (index_range.lower > 0) // expand range downward
		{
			lte -= h.count_at(index_range.lower);
			if (lte*quantile.den < lte_ratio) break;
			--index_range.lower;
		}
		while (index_range.upper+1 < size) // expand range upward
		{
			gte -= h.count_at(index_range.upper);
			if (gte*quantile.den < gte_ratio) break;
			samples_lower += h.count_at(index_range.upper);
			++index_range.upper;
		}
	}
}