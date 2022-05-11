#pragma once

#include <array>
#include <vector>
#include <type_traits>
#include <limits>


namespace quern
{
	// Defined in slice.hpp
	struct array_slice;
	template<size_t N> struct grid_slice;

	/*
		Base type for data grids of arbitrary type and dimensionality
	*/
	class grid_base
	{
	public:
		using index_t = ptrdiff_t;

		static const index_t REJECT = -1;

		// Sampling policy for out-of-range coordinates
		enum OUT_OF_RANGE_POLICY
		{
			// When attempting to query out-of-range coordinates...
			OOR_UNSAFE = 0, // ...undefined behavior
			OOR_FAIL   = 1, // ...fail
			OOR_CLAMP  = 2, // ...clamp to the available range
			OOR_WRAP   = 3, // ...wrap around to the available range
			// TODO option for substitute value?
		};
	};

	/*
		An N-dimensional grid of values, used in data binning.
	*/
	template<typename Value, size_t Dimensionality>
	class grid : public grid_base
	{
	public:
		static const size_t dimensionality = Dimensionality;
		static const size_t N = Dimensionality;

		// Types
		using value_t = Value;
		using coord_t = std::array<index_t, N>;
		//    index_t

		template<typename T_Frac>
		using coord_frac_t = std::array<T_Frac, N>;

		// Filter type
		using filter_t = grid_slice<N>;

		// STL-style aliases
		using value_type = value_t;
		using key_type   = coord_t;

	private:
		// Implementation
		friend class const_iterator;
		using _store_t = std::vector<value_t>;

	public:
		

		// Default interpolator
		template<typename T_Frac>
		struct interpolator_default
		{
			value_t operator()(const value_t &l, const value_t &r, const T_Frac frac) const
			{
				return l + (r - l) * frac;
			}
		};

		// Signifier for constructing end-iterators
		enum iterator_end_t {iterator_end};

		// Iterator implementation
		struct const_iterator
		{
		protected:
			const grid *_g;
			const value_t  *_i;
			coord_t         _c;

			friend class grid;
			const_iterator(const grid *g, const value_t *i, const coord_t &c)
				: _g(g), _i(i), _c(c) {}

			void _inc()
			{
				++_i;
				auto &size = _g->dimensions();
				for (auto d = dimensionality; d--;)
				{
					if (++_c[d] < size[d]) break;
					_c[d] = 0;
				}
			}
			void _dec()
			{
				--_i;
				auto &size = _g->dimensions();
				for (auto d = dimensionality; d--;)
				{
					if (_c[d]-- != 0) break;
					_c[d] = size[d]-1;
				}
			}

		public:
			const_iterator()                                     : _g(nullptr), _i(nullptr), _c{} {}
			const_iterator(const grid &g)                    : _g(&g), _i( g._store.data()), _c{} {}
			const_iterator(const grid &g, iterator_end_t)    : _g(&g), _i(&g._store.back()+1), _c{}
				{_c[dimensionality-1] = _g->dimensions()[dimensionality-1];}

			// Dereference value
			const value_t &operator* () const    {return *_i;}
			const value_t *operator->() const    {return  _i;}

			// Get index or coordinate of this bin
			index_t        index () const    {return _i - _g->_store.data();}
			const coord_t &coord () const    {return _c;}

			// Comparison
			bool operator==(const const_iterator &o) const    {return _i == o._i;}
			bool operator!=(const const_iterator &o) const    {return _i != o._i;}
			bool operator< (const const_iterator &o) const    {return _i <  o._i;}
			bool operator<=(const const_iterator &o) const    {return _i <= o._i;}
			bool operator> (const const_iterator &o) const    {return _i >  o._i;}
			bool operator>=(const const_iterator &o) const    {return _i >= o._i;}

			// Arithmetic
			index_t         operator-(const const_iterator &o) const      {return _i - o._i;}
			const_iterator  operator++(int)                       {auto r=*this; _inc(); return r;}
			const_iterator  operator--(int)                       {auto r=*this; _dec(); return r;}
			const_iterator &operator++()                          {_inc(); return *this;}
			const_iterator &operator--()                          {_dec(); return *this;}
		};

		struct iterator : public const_iterator
		{
		protected:
			friend class grid;
			iterator(grid *g, const value_t *i, const coord_t &c)   : const_iterator(g,i,c) {}

		public:
			iterator()                               : const_iterator()                {}
			iterator(grid &g)                    : const_iterator(g)               {}
			iterator(grid &g, iterator_end_t)    : const_iterator(g, iterator_end) {}

			// Access mutable key
			value_type &operator* () const    {return const_cast<value_type&>(*const_iterator::_i);}
			value_type *operator->() const    {return const_cast<value_type*>( const_iterator::_i);}

			// Arithmetic
			iterator  operator++(int)                       {auto r=*this; this->_inc(); return r;}
			iterator  operator--(int)                       {auto r=*this; this->_dec(); return r;}
			iterator &operator++()                          {this->_inc(); return *this;}
			iterator &operator--()                          {this->_dec(); return *this;}
		};


	public:
		/*
			This default constructor creates a grid with zero elements.
				A reformat will be necessary to get use out if it.
		*/
		grid() : _dims{} {}

		/*
			Set up a uniform grid based on dimensions and initial value.
		*/
		grid(const coord_t &dimensions, const value_t &fill = value_t{})
			: _dims(dimensions) {_store.resize(TotalItems(dimensions), fill);}

		/*
			Clear the grid to the given fill-value.
		*/
		void clear(const value_t &fill = value_t{})
		{
			for (auto &i : _store) i = fill;
		}

		/*
			Reformat the Grid to a new size, erasing all data.
		*/
		void reformat(const coord_t &dimensions, const value_t &fill = value_t{})
		{
			_dims = dimensions;
			_store.clear();
			_store.resize(TotalItems(dimensions), fill);
		}

		/*
			Get the number of items in a grid of the given size.
		*/
		static index_t TotalItems(const coord_t &dimensions)
		{
			index_t n = 1;
			for (size_t d = 0; d < dimensionality; ++d)
			{
				n *= dimensions[d];
				if (n <= 0) return 0;
			}
			return n;
		}

		

		/*
			Access the dimensions
		*/
		size_t           total_size() const    {return _store.size();}
		const coord_t   &dimensions() const    {return _dims;}
		
		/*
			Iterators.
		*/
		const_iterator begin() const    {return const_iterator(*this);}
		iterator       begin()          {return iterator      (*this);}
		const_iterator end  () const    {return const_iterator(*this, iterator_end);}
		iterator       end  ()          {return iterator      (*this, iterator_end);}


		/*
			Get an iterator pointing to the given coordinate or index.
		*/
		const_iterator to      (const coord_t &coord) const    {return const_iterator(this, _store.data()+coord_to_index(coord, _store.size()), coord);}
		/* */ iterator to      (const coord_t &coord)          {return       iterator(this, _store.data()+coord_to_index(coord, _store.size()), coord);}
		const_iterator to_index(const index_t  index) const    {return contains_index(index) ? const_iterator(this, _store.data()+index, index_to_coord(index)) : end();}
		/* */ iterator to_index(const index_t  index)          {return contains_index(index) ?       iterator(this, _store.data()+index, index_to_coord(index)) : end();}

		/*
			Access elements at the given coordinate or index.
		*/
		const value_t &at      (const coord_t &coord, const value_t &out_of_range_value) const    {auto i=coord_to_index(coord); return (i<0) ? out_of_range_value : _store[i];}
		value_t       &at      (const coord_t &coord,       value_t &out_of_range_value)          {auto i=coord_to_index(coord); return (i<0) ? out_of_range_value : _store[i];}
		const value_t &at_index(const index_t  index, const value_t &out_of_range_value) const    {return contains_index(index) ? _store[index] : out_of_range_value;}
		value_t       &at_index(const index_t  index,       value_t &out_of_range_value)          {return contains_index(index) ? _store[index] : out_of_range_value;}

		/*
			Fast, unsafe element access.
		*/
		const value_t &at_unsafe      (const coord_t &coord) const    {return _store[coord_to_index_unsafe(coord)];}
		value_t       &at_unsafe      (const coord_t &coord)          {return _store[coord_to_index_unsafe(coord)];}
		const value_t &at_index_unsafe(const index_t  index) const    {return _store[index];}
		value_t       &at_index_unsafe(const index_t  index)          {return _store[index];}


		/*
			Sample values from the grid...
				sample_index(2) : grab sample at an index
				sample      (2) : grab sample at a coordinate
				sample_clamp(1) : grab sample at closest existing coordinate
				sample_wrap (1) : grab sample, wrapping coordinate if out of range

			sample_clamp and sample_wrap are unsafe on zero-element grids.
		*/
		template<OUT_OF_RANGE_POLICY T_OOR = OOR_FAIL>
		value_t sample      (const coord_t &coord, const value_t out_of_range_value) const    {auto i=coord_to_index<T_OOR>(coord); return (i<0) ? out_of_range_value : _store[i];}
		value_t sample_clamp(const coord_t &coord)                                   const    {return _store[coord_to_index_clamp(coord)];}
		value_t sample_wrap (const coord_t &coord)                                   const    {return _store[coord_to_index_wrap (coord)];}
		value_t sample_index(const index_t  index, const value_t out_of_range_value) const    {return contains_index(index) ? _store[index] : out_of_range_value;}

		/*
			Sample values with a fractional coordinate.
		*/
		template<
			OUT_OF_RANGE_POLICY T_OOR          = OOR_FAIL,
			typename            T_Frac         = float,
			typename            T_Interpolator = interpolator_default<T_Frac>>
		value_t sample(
			coord_frac_t<T_Frac>  coord_frac,
			const value_t         out_of_range_value,
			const T_Interpolator &interpolator = T_Interpolator()) const
		{
			// Calculate floor, ceiling and fractional coordinate
			coord_t cl, ch;
			for (size_t i = 0; i < N; ++i)
			{
				T_Frac floor = std::floor(coord_frac[i]);
				cl[i] = floor;
				ch[i] = std::ceil (coord_frac[i]);
				coord_frac[i] -= floor;
			}

			// Range policy...
			if (T_OOR == OOR_CLAMP || T_OOR == OOR_WRAP)
				{_coord_fix<T_OOR>(cl); _coord_fix<T_OOR>(ch);}
			else if (T_OOR != OOR_UNSAFE && (!contains_coord(cl) || !contains_coord(ch)))
				return out_of_range_value;

			// Preprocess coordinates for sampling subroutine
			index_t m = 1;
			for (size_t i = N-1; i--;) {m *= _dims[i+1]; cl[i] *= m; ch[i] *= m;}

			return _sample_sub<0>(cl, ch, coord_frac, 0, interpolator);
		}
		


		/*
			Convert between coordinates and indices.
				Out-of-range coordinates will yield <on_fail> defaulting to -1.
		*/
		template<OUT_OF_RANGE_POLICY T_OOR = OOR_FAIL>
		index_t coord_to_index      (const coord_t &coord, const index_t on_fail = REJECT) const
		{
			index_t i = 0;
			for (size_t d = 0; d < dimensionality; ++d)
			{
				index_t c = coord[d];
				if      (T_OOR == OOR_WRAP )  c %= _dims[d];
				else if (T_OOR == OOR_CLAMP)  c = std::min<index_t>(std::max<index_t>(c, 0), _dims[d]-1);
				else if (T_OOR != OOR_UNSAFE) if (c < 0 || c >= _dims[d]) return on_fail;
				i = i * _dims[d] + c;
			}
			return i;
		}
		index_t coord_to_index_clamp (const coord_t &coord) const    {return coord_to_index<OOR_CLAMP >(coord);}
		index_t coord_to_index_wrap  (const coord_t &coord) const    {return coord_to_index<OOR_WRAP  >(coord);}
		index_t coord_to_index_unsafe(const coord_t &coord) const    {return coord_to_index<OOR_UNSAFE>(coord);}

		coord_t index_to_coord(index_t index) const
		{
			coord_t c;
			if (contains_index(index))
			{
				for (size_t d = dimensionality; d-- > 0;)
				{
					c[d] = index % _dims[d];
					index /= _dims[d];
				}
			}
			else for (auto &cv : c) cv = REJECT;
			return c;
		}

		/*
			Check if a given index or coordinate is in range.
		*/
		bool contains_index(index_t        index) const
		{
			return index >= 0 && index < index_t(total_size());
		}
		bool contains_coord(const coord_t &coord) const
		{
			for (size_t d = 0; d < dimensionality; ++d)
				if (coord[d] < 0 || coord[d] >= _dims[d]) return false;
			return true;
		}


	private:
		template<OUT_OF_RANGE_POLICY T_OOR>
		void _coord_fix(coord_t &c) const
		{
			for (size_t d = 0; d < N; ++d)
			{
				if (T_OOR == OOR_CLAMP) c[d] = std::min<index_t>(std::max<index_t>(c[d], 0), _dims[d]);
				if (T_OOR == OOR_WRAP ) c[d] %= _dims[d];
			}
		}

		// Interpolated sampling subroutine...
		template<
			size_t   I,
			typename T_Frac,
			typename T_Interpolator>
		value_t _sample_sub(
			const coord_t              &cl,
			const coord_t              &ch,
			const coord_frac_t<T_Frac> &frac,
			index_t                     index,
			const T_Interpolator       &inter) const
		{
			static constexpr size_t I_Next = std::min(I+1, N-1);
			static constexpr bool   LAST   = (I == N-1);

			index += cl[I];
			const value_t &a = (LAST ? _store[index] : _sample_sub<I_Next>(cl,ch,frac,index,inter));
			if (cl[I] == ch[I]) return a;
			index += ch[I]-cl[I];
			const value_t &b = (LAST ? _store[index] : _sample_sub<I_Next>(cl,ch,frac,index,inter));

			return inter(a, b, frac[I]);
		}


	private:
		// Dimensions
		coord_t  _dims;
		_store_t _store;
	};
}
