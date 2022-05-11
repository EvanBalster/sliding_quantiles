#pragma once

#pragma once

#include "quantile.hpp"
#include "binning.hpp"


/*
	binning rules for univariate histograms.

	Data structures for statistical analysis of univariate and multivariate sets.

		Use "Bins" for large-volume collection of sample data.  (millions of samples)
		Use "Scatter" for small-volume collection of sample data, and to prepare binning.

		We can store outliers as scatter data when binning.
*/

namespace quern
{
	/*
		binning_auto_:  automatically computes binning schemes from data.
	*/
	template<typename Quantile = double>
	struct binning_auto_
	{
	public:
		using quantile_t = Quantile;

		/*
			Configuration:
				bins:            number of bins for continuous values
				quantile_min~Max: quantile range for continuous values
		*/
		size_t     bins = 512;
		quantile_t quantile_min = .005, quantile_max = .995;
		
	public:
		binning_auto_(size_t _bins = 512, quantile_t quantileTrim = quantile_t(.005))    : bins(_bins), quantile_min(quantileTrim), quantile_max(quantile_t(1.0)-quantileTrim) {}

		binning_auto_(size_t _bins, quantile_t _quantile_min, quantile_t _quantile_max)    : bins(_bins), quantile_min(_quantile_min), quantile_max(_quantile_max) {}

	
		/*
			Automatic binning for primitive (1DOF) datapoints
		*/
		template<typename DataSet, typename DataPoint = std::decay_t<decltype(*std::declval<const DataSet&>().begin())> >
		std::enable_if_t<std::is_same<DataPoint, bool>::value, binning_params<DataPoint>>
			binning(const DataSet &data) const
		{
			// binning booleans is trivial.
			return {};
		}

		template<typename DataSet, typename DataPoint = std::decay_t<decltype(*std::declval<const DataSet&>().begin())> >
		std::enable_if_t<dof_is_primitive_discrete<DataPoint> && !std::is_same<DataPoint, bool>::value, binning_params<DataPoint>>
			binning(const DataSet &data) const
		{
			// Ignore quantiles when binning discrete values.
			auto range = find_set_range(data);
			return {range.first, range.second};
		}

		template<typename DataSet, typename DataPoint = std::decay_t<decltype(*std::declval<const DataSet&>().begin())> >
		std::enable_if_t<dof_is_primitive_continuous<DataPoint>, binning_params<DataPoint>>
			binning(const DataSet &data) const
		{
			if (quantile_min <= 0.0 && quantile_max >= 1.0)
			{
				// Bin continuous values across full range (use cheaper algorithm)
				auto range = find_set_range(data);
				return {range.first, range.second, bindex_t(bins)};
			}
			else
			{
				// Bin continuous values by quantile
				if (quantile_min >= quantile_max) throw;
				return {
					find_set_quantile(data, quantile_min),
					find_set_quantile(data, quantile_max),
					bindex_t(bins)};
			}
		}
		
		/*
			Automatic binning for multi-dimensional data
		*/
		template<typename DataSet, typename DataPoint = std::decay_t<decltype(*std::declval<const DataSet&>().begin())> >
		std::enable_if_t<!dof_is_primitive<DataPoint>, binning_params<DataPoint>>
			binning(const DataSet &data) const
		{
			return detail::binning_indices<dof_elems<DataPoint>>::
				template map<binning_params<DataPoint>>(*this, data);
		}
	};
	
	template<typename DataSet, typename Quantile = double>
	auto binning_auto(
		const DataSet &data,
		size_t         bins,
		Quantile       quantile_min,
		Quantile       quantile_max)
		-> decltype(std::declval<binning_auto_<Quantile>>().binning(data))
	{
		binning_auto_<Quantile> rule(bins, quantile_min, quantile_max);
		return rule.binning(data);
	}
	
	template<typename DataSet, typename Quantile = double>
	auto binning_auto(
		const DataSet &data,
		size_t         bins         = 512,
		Quantile       quantileTrim = .005)
		-> decltype(std::declval<binning_auto_<Quantile>>().binning(data))
	{
		binning_auto_<Quantile> rule(bins, quantileTrim);
		return rule.binning(data);
	}
}



