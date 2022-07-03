#pragma once


#include "quantile.hpp"
#include "bin_table.hpp"


namespace quern
{
	/*
		Signifier: erase binned samples when creating a histogram.
	*/
	enum EraseBinnedSamples_t {EraseBinnedSamples};


	/*
		A collection of bins quantifying the number of samples in each bin's range.
			"Count" is flexible; float or signed values are permissible.
			Count may be used to represent a weighted or interpolated value.
	*/
	template<
		typename Sample,
		typename Count = uint32_t,
		typename Binning = binning<Sample> >
	class histogram :
		public bin_table<Sample, Count, Binning>
	{
	public:
		using table_t = bin_table<Sample, Count, Binning>;

		using sample_t       = Sample;
		using count_t        = Count;
		using index_t        = typename table_t::index_t;
		using coord_t        = typename table_t::coord_t;
		using binning_t      = typename table_t::binning_t;
		using params_t       = typename binning_t::params_t;
		using iterator       = typename table_t::iterator;
		using const_iterator = typename table_t::const_iterator;

		using quantile_range_t = quantile_range<sample_t>;

		static_assert(
			std::is_integral<count_t>::value && std::is_unsigned<count_t>::value,
			"Bins count type must be unsigned integer.");

	public:
		/*
			Default constructor.  We won't be able to add samples...
		*/
		explicit histogram()    : bin_table<Sample, Count, Binning>() {}

		/*
			Set up empty bins based on an array of binning rules.
		*/
		histogram(const binning_t &binning)    : table_t(binning, count_t(0)) {}
		histogram(const params_t  &params )    : table_t(params , count_t(0)) {}


		/*
			Add or subtract samples.
				Returns whether the sample was in binning range.
		*/
		void add_at(const index_t   index,  const count_t n = 1) noexcept    {count_t dummy; this->at_index(index, dummy) += n;}
		void sub_at(const index_t   index,  const count_t n = 1) noexcept    {count_t dummy; this->at_index(index, dummy) -= n;}
		void add_at(const coord_t  &coord,  const count_t n = 1) noexcept    {return add_at(this->coord_to_index(coord), n);}
		void sub_at(const coord_t  &coord,  const count_t n = 1) noexcept    {return sub_at(this->coord_to_index(coord), n);}
		void add   (const sample_t &sample, const count_t n = 1) noexcept    {return add_at(this->index_for(sample), n);}
		void sub   (const sample_t &sample, const count_t n = 1) noexcept    {return sub_at(this->index_for(sample), n);}
		

		/*
			Access or increment the count at the given indices.
		*/
		// /* */ count_t &count_at(const index_t  i)          {return this->at_index(i,0);}
		count_t count_at(const index_t  i) const    {return this->at_index(i,0);}
		// /* */ count_t &count_at(const coord_t &c)          {return this->at_coord(c,0);}
		count_t count_at(const coord_t &c) const    {return this->at_coord(c,0);}


		/*
			Calculate the total population by iterating over the histogram.
				Use tracked_histogram for inexpensive access to the total.
		*/
		count_t calc_population() const noexcept    {count_t n=0; for (auto &c:this->grid()) n+=c; return n;}

		
#if 0
	public:
		/*
			Automatically configure binning and populate bins from a dataset.
				
				Pass EraseBinnedSamples to remove samples that have been binned.
		*/
		template<typename DataSet, typename Quantile = double>
		histogram(
			const DataSet &data,
			index_t        binsPerDOF   = 1000,
			Quantile       trimQuantile = Quantile(.001)) :
			table_t(binning_auto(data, binsPerDOF, trimQuantile))
		{
			for (auto &sample : data) add(sample);
		}
		template<typename DataSet, typename Quantile = double>
		histogram(
			EraseBinnedSamples_t,
			DataSet  &data,
			index_t   binsPerDOF   = 1000,
			Quantile  trimQuantile = Quantile(.001)) :
			table_t(binning_auto(data, binsPerDOF, trimQuantile))
		{
			auto i = data.begin(), w = i, e = data.end();
			for (; i != e; ++i) if (!add(*i)) *(w++) = *i;
			data.resize(w-data.begin());
		}
#endif
	};



	/*
		Find a quantile in the given histogram.
	*/
	template<typename QuantileInt, typename Sample, typename Count, typename Binning>
	quantile_range<bindex_t> find_quantile_indexes(
		const histogram<Sample, Count, Binning> &histogram,
		const quantile_fraction<QuantileInt>      quantile)
	{
		static_assert(quern::histogram<Sample,Count,Binning>::dimensionality == 1,
			"find_quantile requires 1D histogram.");

		Count numerator = quantile.num, denominator = quantile.den;

		Count quota = histogram.calc_population() * numerator, leq = histogram.count_at(0)*denominator;
		bindex_t size = histogram.bins(), index = 0;

		while (index+1 < size && leq < quota) leq += histogram.count_at(++index)*denominator;

		quantile_range<bindex_t> result;
		result.lower = index;
		if (leq == quota)
			while (index+1 < size) {if (histogram.count_at(++index)) break;}
		result.upper = index;
		return result;

#if 0
		quantile_range<Sample> result = {histogram.binning().min(), histogram.binning().max()};

		bindex_t i_low = 0, i_high = histogram.bins();
		Count low = histogram.count_at(i_low), high = histogram.count_at(i_high);

		while (i_low < i_high)
		{
			bindex_t advance = i_low;
			while (low < high) low += histogram.count_at(++advance);
			i_low = advance;
			advance = i_high;
			while (high < low) high += histogram.count_at(--advance);
			i_high = advance;
		}
#endif
	}

	template<typename QuantileInt, typename Sample, typename Count, typename Binning>
	quantile_range<Sample> find_quantile(
		const histogram<Sample, Count, Binning> &histogram,
		const quantile_fraction<QuantileInt>      quantile)
	{
		auto indexes = find_quantile_indexes(histogram, quantile);
		auto &rule = histogram.binning();
		return {rule.min({indexes.lower}), rule.max({indexes.upper})};
	}
}
