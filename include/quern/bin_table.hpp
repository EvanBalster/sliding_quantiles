#pragma once

#include "binning.hpp"
#include "grid.hpp"

namespace quern
{
	/*
		A table of values keyed by multidimensional binnable variables.
			uses grid<Value> as the internal data store.
	*/
	template<
		typename Key,
		typename Value,
		typename Binning = binning<Key> >
	class bin_table : public grid_base
	{
	public:
		static_assert(dof_count<Key>, "bin_table key must have at least one degree of freedom.");

		// Binning
		using key_t     = Key;
		using binning_t = Binning;
		using params_t  = typename binning_t::params_t;

		// Base grid type
		using grid_t    = grid<Value, dof_count<Key>>;
		using value_t   = typename grid_t::value_t;
		using index_t   = typename grid_t::index_t;
		using coord_t   = typename grid_t::coord_t;

		// STL style aliases
		using key_type   = key_t;
		using value_type = value_t;

		// Dimensionality
		static const size_t dimensionality = grid_t::dimensionality;

		// Fractional coordinate
		template<typename T_Frac>
		using coord_frac_t = typename grid_t::template coord_frac_t<T_Frac>;


	public:
		// Signifier for constructing end-iterators
		using iterator_end_t = typename grid_t::iterator_end_t;
		static const iterator_end_t iterator_end = grid_t::iterator_end;

		using base_iterator = typename grid_t::const_iterator;

		// Iterator implementation
		struct const_iterator : public grid_t::const_iterator
		{
		protected:
			friend class bin_table;
			using base_t = typename grid_t::const_iterator;

			const bin_table *_b;

		public:
			const_iterator()                                      : base_t(), _b(nullptr) {}
			const_iterator(const bin_table &b, const base_t &i)   : base_t(i), _b(&b) {}
			const_iterator(const bin_table &b)                    : base_t(b._grid), _b(&b) {}
			const_iterator(const bin_table &b, iterator_end_t)    : base_t(b._grid,iterator_end), _b(&b) {}

			// Get key range
			key_type key_min() const    {return _b->_binning.min(this->_c);}
			key_type key_max() const    {return _b->_binning.max(this->_c);}
			key_type key_mid() const    {return _b->_binning.mid(this->_c);}
			key_type key    () const    {return _b->_binning.mid(this->_c);}
		};
		struct iterator : public grid_t::iterator
		{
		protected:
			friend class bin_table;
			using base_t = typename grid_t::iterator;

			bin_table *_b;

		public:
			iterator()                                : base_t(), _b(nullptr) {}
			iterator(bin_table &b, const base_t &i)   : base_t(i), _b(&b) {}
			iterator(bin_table &b)                    : base_t(b._grid), _b(&b) {}
			iterator(bin_table &b, iterator_end_t)    : base_t(b._grid,iterator_end), _b(&b) {}

			// Get key range
			key_type key_min() const    {return _b->_binning.min(this->_c);}
			key_type key_max() const    {return _b->_binning.max(this->_c);}
			key_type key_mid() const    {return _b->_binning.mid(this->_c);}
			key_type key    () const    {return _b->_binning.mid(this->_c);}
		};


	public:
		/*
			Default constructor.  We won't be able to add samples...
		*/
		explicit bin_table() {}

		/*
			Set up empty bins based on an array of binning rules.
		*/
		bin_table(const binning_t &binning, const value_t &fill = value_t{})
			: _grid(binning.grid_size(), fill), _binning(binning) {}

		/*
			Clear all values in the BinMap.
		*/
		void clear(const value_t &fill = value_t{})
		{
			_grid.clear(fill);
		}

		/*
			Reformat the BinMap with a new binning rule, erasing all data.
		*/
		void reformat(const binning_t &binning, const value_t &fill = value_t{})
		{
			_binning = binning;
			_grid.reformat(binning.grid_size(), fill);
		}

		/*
			Access the underlying data grid.
		*/
		const grid_t &grid() const    {return _grid;}


		/*
			Access the binning scheme and total number of bins.
		*/
		coord_t          dimensions() const    {return _grid.dimensions();}
		coord_t          grid_size()  const    {return _grid.dimensions();}
		index_t          bins()       const    {return _grid.total_size();}
		const binning_t &binning()    const    {return _binning;}

		index_t coord_to_index(const coord_t &coord) const    {return _grid.coord_to_index(coord);}
		coord_t index_to_coord(const index_t &index) const    {return _grid.index_to_coord(index);}

		
		/*
			Iterators.
		*/
		const_iterator begin() const    {return const_iterator(*this, _grid.begin());}
		iterator       begin()          {return iterator      (*this, _grid.begin());}
		const_iterator end  () const    {return const_iterator(*this, _grid.end());}
		iterator       end  ()          {return iterator      (*this, _grid.end());}


		/*
			Get the coordinate/index corresponding to a key.
				Could be BIN_REJECT.
		*/
		coord_t        coord_for(const key_t &key) const    {return _binning.coord(key); }
		index_t        index_for(const key_t &key) const    {return coord_to_index(coord_for(key));}

		template<typename R>
		coord_frac_t<R> coord_frac_for(const key_t &key) const    {return _binning.template coord_frac<R>(key);}

		/*
			Access values by key.
		*/
		const_iterator find (const key_t &key) const    {return const_iterator(*this, _grid.to(coord_for(key)));}
		iterator       find (const key_t &key)          {return       iterator(*this, _grid.to(coord_for(key)));}

		/*
			Access elements by key.
		*/
		const value_t &at(const key_t &key,      const value_t &out_of_range_value) const    {return at(coord_for(key), out_of_range_value);}
		value_t       &at(const key_t &key,            value_t &out_of_range_value)          {return at(coord_for(key), out_of_range_value);}
		const value_t &at(const coord_t &coord,  const value_t &out_of_range_value) const    {return _grid.at(coord, out_of_range_value);}
		value_t       &at(const coord_t &coord,        value_t &out_of_range_value)          {return _grid.at(coord, out_of_range_value);}
		const value_t &at_index(const index_t i, const value_t &out_of_range_value) const    {return _grid.at_index(i, out_of_range_value);}
		value_t       &at_index(const index_t i,       value_t &out_of_range_value)          {return _grid.at_index(i, out_of_range_value);}

		/*
			Sample by key.
		*/
		template<
			OUT_OF_RANGE_POLICY T_OOR          = OOR_FAIL,
			typename            T_Frac         = float,
			typename            T_Interpolator = typename grid_t::template interpolator_default<T_Frac>>
		value_t sample(
			const key_t          &key,
			const value_t         out_of_range_value,
			const T_Interpolator &interpolator = T_Interpolator()) const
		{
			auto frac = coord_frac_for<T_Frac>(key);
			return _grid.template sample<T_OOR>(frac, out_of_range_value, interpolator);
		}


		/*
			Get an iterator pointing to the given coordinate or index.
		*/
		const_iterator to_coord(const coord_t &coord) const    {return const_iterator(*this, _grid.to      (coord));}
		/* */ iterator to_coord(const coord_t &coord)          {return       iterator(*this, _grid.to      (coord));}
		const_iterator to_index(const index_t  index) const    {return const_iterator(*this, _grid.to_index(index));}
		/* */ iterator to_index(const index_t  index)          {return       iterator(*this, _grid.to_index(index));}


	private:
		// Binning scheme: can't be changed without resetting store
		binning_t _binning;
		grid_t    _grid;
	};
}
