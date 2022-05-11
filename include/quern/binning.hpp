#pragma once

#pragma once

#include <cmath>
#include <complex>
#include <type_traits>
#include <algorithm>
#include <array>
#include <vector>
#include <stdint.h>

#include "dof.hpp"


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
		binning uses a 64-bit index type.
	*/
	using bin_index_t = ptrdiff_t;
	using bindex_t = bin_index_t;

	/*
		binning coordinates, used for multivariate maps.
	*/
	template<size_t N>                using bin_coord_t      = std::array<bindex_t, N>;
	template<typename Real, size_t N> using bin_coord_frac_t = std::array<Real, N>;

	/*
		This special index is used to indicate that a sample has been rejected from binning.
	*/
	enum BIN_REJECT_t : bindex_t {BIN_REJECT = -1};

	/*
		This structure is used to express the domain of values for a single binning dof.
	*/
	template<typename Real>
	struct bin_domain
	{
		Real min, max;
	};

	/*
		Multivariate domain.
	*/
	template<typename Real, size_t N>
	using grid_domain = std::array<bin_domain<Real>, N>;

	/*
		binning<T> defines a rule for organizing data points into rectangular bins.
			Each degree of freedom in the data introduces a new binning axis.
	*/
	template<class T, class Enable = void> struct binning_params_ {};

	template<class T> using binning_params = typename binning_params_<T>::type;


	template<class T, class Enable = void>
	struct binning
	{
		// (This is an undefined default scheme, used to document the others.)
		static_assert(dof_count<T> >= 0, "Don't know how to bin this type.");

		static const size_t dof = dof_count<T>;
		

		/*
			Each binning scheme defines some instantiation of bin_coords_t as its coord_t.
		*/
		using value_t = T;
		using index_t = bindex_t;
		using coord_t = std::array<bindex_t, 0>;

		/*
			Each binning scheme is defined by some parameters supplied to its constructor.
		*/
		using params_t = binning_params<T>;

		binning();
		binning(const params_t&);

		/*
			Query the total number of bins in this scheme.
		*/
		index_t bins()     const;
		coord_t grid_size() const;

		/*
			Query the domain of the binning scheme.
		*/
		template<typename Real>
		grid_domain<Real, dof> domain() const    {return {};}

		T min() const;
		T max() const;

		/*
			Query the extents of a given bin, and its central value.
		*/
		T min(coord_t  c) const;
		T max(coord_t  c) const;
		T mid(coord_t  c) const;

		/*
			Obtain bin index information for a given value.

			Single-dimension binning schemes may offer additional methods.
		*/
		coord_t  coord (const T&) const; // return grid coordinate; may include BIN_REJECT elements.
		bindex_t accept(const T&) const; // return whether a value can be binned
		bindex_t reject(const T&) const; // return whether a value cannot be binned

		template<typename R>
		bin_coord_frac_t<R, 0> coord_frac(const T&) const;
	};


	// binning for primitive continuous values.
	template<class T>
	struct binning_params_<T, std::enable_if_t<std::is_floating_point<T>::value>>
	{
		using type = binning_params_;

		T        min, max;
		bindex_t bins;

		// Scale resolution
		static type scale(const type &params, bindex_t scale)    {auto p=params; p.bins *= scale; return p;}
	};


	template<class T>
	struct binning<T, std::enable_if_t<std::is_floating_point<T>::value>>
	{
	public:
		static const size_t dof = dof_count<T>;

		using value_t = T;
		using index_t = bindex_t;
		using coord_t = bin_coord_t<1>;
		using params_t = binning_params<T>;

	public:
		// Default constructor: no binning
		binning() : _min(0.0), _max(0.0), _step(1.0), _bins(0) {}

		// Constructor
		binning(const params_t &p) :
			_min(p.min), _max(p.max),
			_step((p.max-p.min)/T(std::max(p.bins, bindex_t(1)))),
			_bins(                std::max(p.bins, bindex_t(1))) {}

		// Get parameters
		params_t params() const    {return {_min, _max, _bins};}

		// Get extents
		template<typename Real>
		grid_domain<Real, 1> domain() const    {return {{_min, _max}};}

		T min()          const    {return _min;}
		T max()          const    {return _max;}
		T step()         const    {return _step;}           // step() is unique to this type
		T min(coord_t c) const    {return _min + _step * c[0];}
		T max(coord_t c) const    {return min(c) + _step;}
		T mid(coord_t c) const    {return min(c) + _step * T(.5);}

		// Grid size
		index_t bins()     const    {return _bins;}
		coord_t grid_size() const    {return {bins()};}

		// binning queries.
		bool     accept(const T v) const    {return v >= _min && v <  _max;}
		bool     reject(const T v) const    {return v <  _min || v >= _max;}
		coord_t  coord (const T v) const    {return {index(v)};}
		bindex_t index (const T v) const
		{
			return reject(v) ? BIN_REJECT : std::min(_vi(v), _bins-1);
		}

		// Real-valued coordinate
		template<typename R>
		bin_coord_frac_t<R, 1> coord_frac(const T v) const    {return {(v-_min)/_step - R(.5)};}


	private:
		T        _min = 0.0, _max = 0.0, _step = 1.0;
		bindex_t _bins = 0;
		
		// subroutines
		index_t _vi(const T v) const    {return index_t((v-_min)/_step);}
	};

	// binning for booleans.
	template<class T>
	struct binning_params_<T, std::enable_if_t<std::is_same<T,bool>::value>>
	{
		using type = binning_params_;

		// Scale resolution (no effect)
		static type scale(const type &params, bindex_t scale)    {return params;}
	};

	template<class T>
	struct binning<T, std::enable_if_t<std::is_same<T,bool>::value>>
	{
	public:
		static const size_t dof = dof_count<T>;

		using value_t = T;
		using index_t = bindex_t;
		using coord_t = bin_coord_t<1>;
		using params_t = binning_params<T>;

	public:
		// Default constructor: single bin at zero
		binning() {}

		// Constructor
		binning(const params_t &p) {}

		// Get parameters
		params_t params() const    {return {};}

		// Extents
		template<typename Real>
		grid_domain<Real, 1> domain() const    {return {{-0.5, 1.5}};}

		T min()          const    {return false;}
		T max()          const    {return true;}
		T min(coord_t c) const    {return c[0] > 0;}
		T max(coord_t c) const    {return c[0] > 0;}
		T mid(coord_t c) const    {return c[0] > 0;}

		// Grid size
		index_t bins()     const    {return 2;}
		coord_t grid_size() const    {return {2};}

		// binning queries.
		bool    accept(const T v) const    {return true;}
		bool    reject(const T v) const    {return false;}
		coord_t coord (const T v) const    {return {index(v)};}
		index_t index (const T v) const    {return v ? 1 : 0;}

		// Real-valued coordinate
		template<typename R>
		bin_coord_frac_t<R, 1> coord_frac(const T v) const    {return {v ? R(1) : R(0)};}
	};


	/*
		Binning for primitive discrete values (integers and enums).
			Allowed values should be consecutive.
	*/
	template<class T>
	struct binning_params_<T, std::enable_if_t<dof_is_primitive_discrete<T> && !std::is_same<T,bool>::value>>
	{
		using type = binning_params_;

		// Inclusive range of consecutive values.
		T min, max;

		// Scale resolution (no effect)
		static type scale(const type &params, bindex_t scale)    {return params;}
	};

	template<class T>
	struct binning<T, std::enable_if_t<dof_is_primitive_discrete<T> && !std::is_same<T,bool>::value>>
	{
	public:
		static const size_t dof = dof_count<T>;

		using value_t = T;
		using index_t = bindex_t;
		using coord_t = bin_coord_t<1>;
		using params_t = binning_params<T>;

	public:
		// Default constructor: single bin at zero
		binning() : _min{}, _max{} {}

		// Constructor
		binning(const params_t &p) : _min(p.min), _max(p.max) {}

		// Get parameters
		params_t params() const    {return {T(_min), T(_max)};}

		// Extents
		template<typename Real>
		grid_domain<Real, 1> domain() const    {return {{_min-Real(0.5), _max+Real(0.5)}};}

		T min()          const    {return T(_min);}
		T max()          const    {return T(_max);}
		T min(coord_t c) const    {return T(_min+c[0]);}
		T max(coord_t c) const    {return T(_min+c[0]);}
		T mid(coord_t c) const    {return T(_min+c[0]);}

		// Grid size
		index_t bins()     const    {return (_max-_min)+1;}
		coord_t grid_size() const    {return {bins()};}

		// binning queries.
		bool    accept(const T v) const    {return index_t(v) >= _min && index_t(v) <= _max;}
		bool    reject(const T v) const    {return index_t(v) <  _min || index_t(v) >  _max;}
		coord_t coord (const T v) const    {return {index(v)};}
		index_t index (const T v) const
		{
			return reject(v) ? BIN_REJECT : _vi(v);
		}

		// Real-valued coordinate
		template<typename R>
		bin_coord_frac_t<R, 1> coord_frac(const T v) const    {return {v-_min};}

	private:
		bindex_t _min, _max;
		
		// subroutines
		index_t _vi(const T v) const    {return index_t(v)-_min;}
	};
}



