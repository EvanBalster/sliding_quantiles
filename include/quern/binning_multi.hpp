#pragma once

#pragma once

#include <cmath>
#include <complex>
#include <tuple>
#include <functional>
#include <type_traits>
#include <algorithm>
#include <array>
#include <vector>
#include <stdint.h>

#include "binning.hpp"


/*
	Binning rules for multivariate histograms.

	Data structures for statistical analysis of univariate and multivariate sets.

		Use "Bins" for large-volume collection of sample data.  (millions of samples)
		Use "Scatter" for small-volume collection of sample data, and to prepare binning.

		We can store outliers as scatter data when binning.
*/

namespace quern
{
	namespace detail
	{
		template<size_t... I>
		struct binning_indices_
		{
			using coord_t = std::array<bindex_t, sizeof...(I)>;

			template<class Real, size_t dof, class B>
			static grid_domain<Real, dof> domain(const B &b)                      {return {std::get<I>(b).template domain<Real>()[0] ...};}

			template<class V, class B> static V min(const B &b)                      {return V(std::get<I>(b).min() ...);}
			template<class V, class B> static V max(const B &b)                      {return V(std::get<I>(b).max() ...);}
			template<class V, class B> static V min(const B &b, const coord_t &c)    {return V(std::get<I>(b).min(c[I]) ...);}
			template<class V, class B> static V max(const B &b, const coord_t &c)    {return V(std::get<I>(b).max(c[I]) ...);}
			template<class V, class B> static V mid(const B &b, const coord_t &c)    {return V(std::get<I>(b).mid(c[I]) ...);}
			
			template<class Dest, class Src>
			static constexpr Dest make(const Src &params)
				{return Dest(std::get<I>(params) ...);}

			template<class Dest, class Src>
			static constexpr Dest params(const Src &binnings)
				{return Dest(std::get<I>(binnings).params() ...);}

			template<class Key, class Params> // Ugh
			static constexpr Params scale(const Params &params, bindex_t scale)
				{return Params(binning_params_<typename std::tuple_element<I,Key>::type>::scale(std::get<I>(params), scale) ...);}

			template<class Dest, class Func, class Src>
			static constexpr Dest map(const Func &func, const Src &params)
				{return {func(std::get<I>(params)) ...};}
		};
		
		template<typename T>                  struct binning_indices_make;
		template<typename Int_t, size_t... I> struct binning_indices_make<std::integer_sequence<Int_t, I...>> {using type = binning_indices_<I...>;};
		
		template<size_t N> using binning_indices = typename binning_indices_make<std::make_index_sequence<N>>::type;
		
		
		template<typename T>
		struct BinningTuples;
		
		template<class... T_Elems>
		struct BinningTuples<std::tuple<T_Elems...>>
		{
			using indices = binning_indices<std::tuple_size<std::tuple<T_Elems...>>::value>;
			
			using TupleValue = std::tuple<T_Elems...>;
			
			using SubBinning = std::tuple<binning<T_Elems>        ...>;
			using SubParams  = std::tuple<binning_params<T_Elems> ...>;
			
			static SubBinning MakeBinning(const SubParams &p)    {return indices::template make<SubBinning>(p);}

			static SubParams GetParams(const SubBinning &p)    {return indices::template params<SubParams>(p);}
		};
	}
	
	/*
		binning for non-primitive values with one or more degrees of freedom
	*/
	template<class T>
	struct binning_params_<T, std::enable_if_t<!dof_is_primitive<T> && (dof_count<T> > 0)>>
	{
		using _tuples  = detail::BinningTuples<typename DOF_Info<T>::tuple_t>;
		using _indices = typename _tuples::indices;


		using type = typename _tuples::SubParams;

		static type scale(const type &params, bindex_t scale)
			{return _indices::template scale<T, type>(params, scale);}
	};

	template<class T>
	struct binning<T, std::enable_if_t<!dof_is_primitive<T> && (dof_count<T> > 0)>>
	{
	private:
		// Implementation details
		using _tuples = detail::BinningTuples<typename dof_info<T>::tuple_t>;
	
	public:
		static const size_t dof = dof_count<T>;

		using value_t = T;
		using index_t = bindex_t;
		using coord_t = bin_coord_t<dof_count<T>>;
		using params_t = binning_params<T>;
		
	public:
		// Default constructor (default binning per dof)
		binning() {}

		// Constructor
		binning(const params_t &p) : _sub(_tuples::MakeBinning(p)) {}

		// Get parameters
		params_t params() const    {return _tuples::GetParams(_sub);}
		
		// Extents
		template<typename Real>
		grid_domain<Real, dof> domain() const    {return _tuples::indices::template domain<Real, dof>(_sub);}

		value_t min()          const    {return _tuples::indices::template min<value_t>(_sub);}
		value_t max()          const    {return _tuples::indices::template max<value_t>(_sub);}
		value_t min(coord_t c) const    {return _tuples::indices::template min<value_t>(_sub,c);}
		value_t max(coord_t c) const    {return _tuples::indices::template min<value_t>(_sub,c);}
		value_t mid(coord_t c) const    {return _tuples::indices::template mid<value_t>(_sub,c);}
		
		// Grid size
		index_t bins()     const    {return _bins();}
		coord_t grid_size() const    {coord_t c; _grid_size(c); return c;}

		// binning queries.
		bool    accept(const value_t &v) const    {return _accept(v);}
		bool    reject(const value_t &v) const    {return _reject(v);}
		coord_t coord (const value_t &v) const
		{
			coord_t c; _coord(c,v);
			return c;
		}

		// Real-valued coordinate
		template<typename R>
		bin_coord_frac_t<R, dof> coord_frac(const value_t &v) const
		{
			bin_coord_frac_t<R, dof> c; _coord_frac(c, v);
			return c;
		}
		
	private:
		typename _tuples::SubBinning _sub;

	public:
		// Access binning rules for each element
		template<size_t I> auto element()       -> decltype(std::get<I>(_sub))    {return std::get<I>(_sub);}
		template<size_t I> auto element() const -> decltype(std::get<I>(_sub))    {return std::get<I>(_sub);}

	private:
		// Implementation
		static const size_t _n_elems = dof_info<T>::elems - 1;
		
	#define ELEMENT_SUBROUTINE(RETURN_T, DECL, EXPR, EXPR_ZERO) \
		template<size_t I=_n_elems> std::enable_if_t<(I!=0), RETURN_T> DECL {EXPR;} \
		template<size_t I=_n_elems> std::enable_if_t<(I==0), RETURN_T> DECL {EXPR_ZERO;}
	#define ELEMENT_SUBROUTINE_EX(TEMPLATE, RETURN_T, DECL, EXPR, EXPR_ZERO) \
		template<size_t I=_n_elems, TEMPLATE> std::enable_if_t<(I!=0), RETURN_T> DECL {EXPR;} \
		template<size_t I=_n_elems, TEMPLATE> std::enable_if_t<(I==0), RETURN_T> DECL {EXPR_ZERO;}
		
		ELEMENT_SUBROUTINE(void, _grid_size(coord_t &c) const,
			c[I] = element<I>().bins(); _grid_size<I-1>(c),
			c[I] = element<I>().bins())
		ELEMENT_SUBROUTINE(index_t, _bins() const,
			return element<I>().bins()*_bins<I-1>(),
			return element<I>().bins())
		ELEMENT_SUBROUTINE(bool, _accept(const value_t &v) const,
			return element<I>().accept(dof_elem<I>(v)) && _accept<I-1>(v),
			return element<I>().accept(dof_elem<I>(v)))
		ELEMENT_SUBROUTINE(bool, _reject(const value_t &v) const,
			return element<I>().reject(dof_elem<I>(v)) || _reject<I-1>(v),
			return element<I>().reject(dof_elem<I>(v)))
		ELEMENT_SUBROUTINE(void, _coord(coord_t &c, const value_t &v) const,
			c[I] = element<I>().index(dof_elem<I>(v)); _coord<I-1>(c,v),
			c[I] = element<I>().index(dof_elem<I>(v)))
		ELEMENT_SUBROUTINE_EX(typename R, void, _coord_frac(bin_coord_frac_t<R, dof> &c, const value_t &v) const,
			c[I] = element<I>().template coord_frac<R>(dof_elem<I>(v))[0]; _coord_frac<I-1>(c,v),
			c[I] = element<I>().template coord_frac<R>(dof_elem<I>(v))[0])
	
	#undef ELEMENT_SUBROUTINE
	};
}



