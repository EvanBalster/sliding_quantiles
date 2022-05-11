#pragma once

#include <algorithm>
#include <limits>


namespace quern
{
	/*
		A filter that selects a subset of array elements by one of two methods:
			- A range with minimum and maximum, or
			- A bitmask selecting up to (8*sizeof(ptrdiff_t)-1) items
	*/
	struct array_slice
	{
		enum : ptrdiff_t    {N_MAX = std::numeric_limits<ptrdiff_t>::max()};

		union
		{
			ptrdiff_t lo;
			ptrdiff_t mask;
		};
		ptrdiff_t hi;

		// Return a number with the lowest N bits enabled.
		static constexpr ptrdiff_t _masklet(ptrdiff_t N)
			{return (1 << std::min<ptrdiff_t>(std::max<ptrdiff_t>(N, 0), sizeof(ptrdiff_t)*8-1))-1;}

		// Get the number of elements selected by this filter.
		size_t count() const
		{
			if (is_range()) return std::max<ptrdiff_t>(hi-lo, 0);
			size_t c = 0;
			for (auto m=mask; m; m >>= 1) c += (m&1);
			return c;
		}
		size_t count(size_t array_size) const
		{
			if (is_range()) return Range(lo, std::min<ptrdiff_t>(hi, array_size)).count();
			else            return Mask(mask & _masklet(array_size)).count();
		}

		// Classify this filter
		bool is_range() const    {return hi!=0;}
		bool is_mask () const    {return hi==0;}
		bool is_all  () const    {return lo==0 && hi == N_MAX;}

		bool is_all(ptrdiff_t size) const    {return lo==0 && hi >= size;}

		// Get this filter's mask, or the nearest possible representation for non-mask filters
		ptrdiff_t toMask() const
		{
			if (is_range()) return _masklet(hi) & (~_masklet(lo));
			else return mask;
		}

		// Accept this sample?
		bool accept(ptrdiff_t val) const
		{
			return is_range()
				? (lo <= val && val < hi)
				: (val >= 0 && ((ptrdiff_t(1)<<val) & mask));
		}

		// Set a contiguous, demi-exclusive range to accept
		static array_slice Mask (size_t mask)               {return array_slice{{ptrdiff_t(mask)},0};}
		static array_slice Range(size_t min, size_t max)    {return (max > min) ? array_slice{{ptrdiff_t(min)},ptrdiff_t(max)} : array_slice{{0},0};}
		static array_slice Value(size_t value)              {return array_slice{{ptrdiff_t(value)},ptrdiff_t(value)+1};}
		static array_slice True ()                          {return {{2},0};}
		static array_slice False()                          {return {{1},0};}
		static array_slice All  ()                          {return {{0},N_MAX};}
		static array_slice None ()                          {return {{0},0};}

		// Get the intersection of two filters
		array_slice operator&(const array_slice &o) const
		{
			// Mask intersection
			if (this->is_mask())
			{
				if (o.is_mask()) return Mask(mask & o.mask);
				else             return Mask(mask & o.toMask());
			}
			else if (o.is_mask()) return Mask(toMask() & o.mask);

			// Range intersection
			if (lo > o.hi || o.lo > hi) return None();
			return Range(std::max(lo,o.lo), std::min(hi,o.hi)); 
		}
		array_slice &operator&=(const array_slice &o)
		{
			return *this = (*this & o);
		}

		// Test an array of coordinates against an array of masks
		static bool Accept(const ptrdiff_t *coord, const array_slice *masks, const size_t N)
		{
			for (size_t i = 0; i < N; ++i)
				if (!masks[i].accept(coord[i])) return false;
			return true;
		}
	};

	template<size_t N>
	struct grid_slice
	{
		using coord_t = std::array<ptrdiff_t, N>;
		using array_slices_t = std::array<array_slice, N>;

		array_slices_t filter;

		grid_slice(std::vector<array_slice> list)
		{
			size_t i = 0;
			for (auto &item : list) {filter[i++]=item; if (i>=N) break;}
			while (i < N) filter[i++] = array_slice::All();
		}
		grid_slice(std::initializer_list<array_slice> list)
		{
			size_t i = 0;
			for (auto &item : list) {filter[i++]=item; if (i>=N) break;}
			while (i < N) filter[i++] = array_slice::All();
		}

		// Count the number of bins selected by this filter.
		size_t count(const coord_t &grid_size) const
		{
			size_t total = 1;
			for (size_t i = 0; i < N; ++i)
				total *= (grid_size[i] > 0) ? filter[i].count(grid_size[i]) : size_t(0);
			return total;
		}

		const array_slice &operator[](size_t i) const    {return filter[i];}
		/* */ array_slice &operator[](size_t i)          {return filter[i];}

		bool accept(const ptrdiff_t *coord) const    {return array_slice::Accept(coord, &filter[0], N);}
		bool accept(const coord_t   &coord) const    {return accept(&coord[0]);}

		// Produce a bitmask indicating which filters specify All
		unsigned mask_all() const
		{
			unsigned r = 0;
			for (size_t i = 0; i < N; ++i) if (filter[i].is_all()) r |= (1u << i);
			return r;
		}
		unsigned mask_all(const coord_t &size) const
		{
			unsigned r = 0;
			for (size_t i = 0; i < N; ++i) if (filter[i].is_all(size[i])) r |= (1u << i);
			return r;
		}

		// Iterate over all passing elements in a grid
		template<typename Op>
		void for_each(const coord_t &grid_size, const Op &op) const
		{
			PG_TRACE;

			coord_t coord = {0};
			for_each_sub(grid_size, op, coord, 0);
		}

		// Sum all passing elements in a datagrid
		/*template<typename Value>
		Value sum(const grid<Value, N> &grid) const
		{
			Value sum = 0;
			for_each(grid.dimensions(), [&](const coord_t &coord, ptrdiff_t index) {sum += grid.at_index(index,0);});
			return sum;
		}*/

	private:
		template<typename Op, size_t I = 0>
		void for_each_sub(const coord_t &grid_size, const Op &op, coord_t &coord, ptrdiff_t superIndex) const
		{
			ptrdiff_t baseIndex = superIndex * grid_size[I];

			auto &f = filter[I];
			if (f.hi)
			{
				// Iterate over filter range, without passing grid_size
				auto hi = std::min(f.hi, grid_size[I]);
				for (ptrdiff_t i = f.lo; i < hi; ++i)
				{
					coord[I] = i;
					if (I >= N-1) op((const coord_t&) coord, baseIndex+i);
					else          for_each_sub<Op, (I+1) % N>(grid_size, op, coord, baseIndex+i);
				}
			}
			else
			{
				// Iterate over filter bitmask, without passing grid_size
				auto mask = f.mask & ((1 << std::min<ptrdiff_t>(grid_size[I], 8*sizeof(ptrdiff_t)-1)) - 1);
				for (ptrdiff_t i = 0; mask; mask >>= 1, ++i) if (mask & 1)
				{
					coord[I] = i;
					if (I >= N-1) op((const coord_t&) coord, baseIndex+i);
					else          for_each_sub<Op, (I+1) % N>(grid_size, op, coord, baseIndex+i);
				}
			}
		}
	}
}
