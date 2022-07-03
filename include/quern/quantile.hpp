#pragma once

#include <exception>
#include <utility>
#include <type_traits>
#include <vector>
#include <queue>



namespace quern
{
	/*
		Represents a fraction defining a quantile.
	*/
	template<typename T_Int = ptrdiff_t>
	struct quantile_fraction
	{
		using int_t = T_Int;

		constexpr quantile_fraction(int_t numerator, int_t denominator = 1)
			:
			num(numerator), den(denominator) {}

		union
		{
			struct {int_t num,       den;};
			struct {int_t numerator, denominator;};
		};

		constexpr quantile_fraction operator/(const quantile_fraction &other) const
			{return quantile_fraction(num*other.den, den*other.num);}

		bool operator< (const quantile_fraction &o) const noexcept    {return num*o.den <  o.num*den;}
		bool operator<=(const quantile_fraction &o) const noexcept    {return num*o.den <= o.num*den;}
		bool operator> (const quantile_fraction &o) const noexcept    {return num*o.den >  o.num*den;}
		bool operator>=(const quantile_fraction &o) const noexcept    {return num*o.den >= o.num*den;}
		bool operator==(const quantile_fraction &o) const noexcept    {return num*o.den == o.num*den;}
		bool operator!=(const quantile_fraction &o) const noexcept    {return num*o.den != o.num*den;}

		operator      float () const noexcept    {return num / float(den);}
		operator      double() const noexcept    {return num / double(den);}
		operator long double() const noexcept    {return num / (long double)(den);}
	};

	template<typename T_Int, typename T_IntQuo>
	constexpr auto operator/(T_Int num, const quantile_fraction<T_IntQuo> &den)
			{return quantile_fraction<decltype(T_Int(1)/T_IntQuo(1))>(num*den.den, den.num);}


	namespace literals
	{
		constexpr quantile_fraction<ptrdiff_t> operator""_quo(unsigned long long int num)
		{
			return quantile_fraction<ptrdiff_t>(num,1);
		}
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
		Algorithm for evaluating the range of a non-empty dataset.
	*/
	template<class DataSet, class Value = std::decay_t<decltype(*std::declval<DataSet>().begin())>>
	std::pair<Value, Value> find_set_range(const DataSet &data)
	{
		auto i = data.begin(), e = data.end();
		Value min = *i, max = *i;
		++i;
		while (i != e)
		{
			auto &v = *i;
			if (v < min) min = v;
			if (v > max) max = v;
			++i;
		}
		return std::make_pair(min, max);
	}
	
	/*
		Algorithm for evaluating the quantile function for a non-empty dataset.
			Continuous-valued functions will be interpolated.
	*/
	template<typename Quantile, class DataSet, class Value = std::decay_t<decltype(*std::declval<DataSet>().begin())>>
	Value find_set_quantile(const DataSet &data, const Quantile quantile)
	{
		std::priority_queue<Value, std::vector<Value>, std::greater<Value> > hi;
		std::priority_queue<Value, std::vector<Value>, std::less   <Value> > lo;

		Quantile nLo = 0.0, nTotal = 0.0;
		for (auto i = data.begin(), e = data.end(); i != e; ++i)
		{
			lo.push(*i);
			if (nLo > ++nTotal * quantile) {hi.push(lo.top()); lo.pop();}
			else ++nLo;
		}

		if (std::is_floating_point<Value>::value)
		{
			Quantile mix = nLo - (nTotal*quantile); // nLo-1 < nTotal*quantile <= nLo
			Value vLo = lo.top(), vHi = hi.top();
			return vHi + (vLo - vHi) * mix;
		}
		else return lo.top();
	}
}
