#pragma once


#include <type_traits>
#include <complex>
#include <tuple>


/*
	Compile time degrees-of-freedom analysis.
		Allows supported types to be inspected for individual degrees of freedom.
		Single degrees may be accessed by index.
 
	Supported primitives:
		floating-point values (continuous)
		enumerations (discrete)
 		booleans (discrete)
 
	Supported aggregates:
 		std::complex
		std::tuple
*/


namespace quern
{
	using dof_index_t = size_t;


	// Signifier type used to signal improper DOF
	struct dof_out_of_range {};

	/*
		The dof_count template inspects a type for continuous and discrete degrees of freedom.
			Floating-point values have one continuous degree.
			Booleans and enums have one discrete degree.
			Complex values have twice the degrees of the underlying type.
			Tuples have the combined dimensionality of all elements.
	 
		Each valid dof_info<T> provides the following:
		
		count                -- number of DOF
		elems                -- number of child elements in an aggregate value (0 for primitives)
		primitive            -- whether this type is a primitive DOF value..
		primitive_continuous -- whether this type is a primitive continuous DOF value.
		primitive_discrete   -- whether this type is a primitive discrete DOF value.
 
		value_t              -- decay type of T
		tuple_t           *  -- a tuple comprising the elements of T, if it is an aggregate.
		dof_tuple_t          -- a tuple comprising the primitive degrees of T
		dof_t <N>            -- the type of the Ith degree of freedom in T, if applicable.
		elem_t<N>         *  -- the type of the Ith element in T, if applicable.
 
		dof <N>(T&)          -- access the Nth degree of freedom in T.
		elem<N>(T&)       *  -- access the Nth element of T.
	 
		* These members are only available for supported aggregates (complex & tuple).
	*/
	template<class T, class = void>
	struct dof_info; /* undefined */
	
	
	/*
		Shorthand methods.
			dof_type<N, Value> -- the type of Value's Nth degree of freedom.
			dof_get <N>(value) -- access value's Nth degree of freedom.
	*/
	template<size_t N, class Value> using dof_type     = typename dof_info<Value>::template dof_t<N>;
	template<size_t N, class Value> using dof_elemtype = typename dof_info<Value>::template elem_t<N>;
	
	template<size_t N, class Value> static dof_type    <N, Value> dof_get (Value &v)    {return dof_info<Value>::template dof <N>(v);}
	template<size_t N, class Value> static dof_elemtype<N, Value> dof_elem(Value &v)    {return dof_info<Value>::template elem<N>(v);}
	
	template<class T> static constexpr size_t dof_count                   = dof_info<T>::count;
	template<class T> static constexpr size_t dof_elems                   = dof_info<T>::elems;
	template<class T> static constexpr bool   dof_is_primitive_continuous = dof_info<T>::primitive_continuous;
	template<class T> static constexpr bool   dof_is_primitive_discrete   = dof_info<T>::primitive_discrete;
	template<class T> static constexpr bool   dof_is_primitive            = dof_info<T>::primitive;
	
	
	/*
		The remainder of this header is template implementation.
	*/
	
	template<class T> // Floating-point values
	struct dof_info<T,
		std::enable_if_t<std::is_floating_point<T>::value
		&& !std::is_const<T>::value
		&& !std::is_volatile<T>::value>>
	{
		constexpr static size_t
			count                = 1,
			elems                = 0;
		constexpr static bool
			primitive            = true,
			primitive_continuous = true,
			primitive_discrete   = false;
		
		using value_t            = T;
		using dof_tuple_t        = std::tuple<T>;
		template<size_t N>
		using dof_t              = std::enable_if_t<(N < count), T>;
		
		template<size_t N>
		static dof_t<N> &dof(T &v)    {return v;}
	};

	template<class T> // Enumerations and booleans
	struct dof_info<T,
		std::enable_if_t<(std::is_enum<T>::value || std::is_same<bool, std::remove_cv_t<T>>::value)
			&& !std::is_const<T>::value
			&& !std::is_volatile<T>::value>>
	{
		constexpr static size_t
			count                = 1,
			elems                = 0;
		constexpr static bool
			primitive            = true,
			primitive_continuous = false,
			primitive_discrete   = true;
		
		using value_t            = T;
		using dof_tuple_t        = std::tuple<T>;
		template<size_t N>
		using dof_t              = std::enable_if_t<(N < count), T>;
		
		template<size_t N>
		static dof_t<N> &dof(T &v)    {return v;}
	};

	template<class R> // Complex numbers
	struct dof_info<std::complex<R>>
	{
		constexpr static size_t
			count                = 2*dof_info<R>::count,
			elems                = 2;
		constexpr static bool
			primitive            = false,
			primitive_continuous = false,
			primitive_discrete   = false;
		
		using value_t            = std::complex<R>;
		using tuple_t            = std::tuple<R,R>;
		using dof_tuple_t        = std::tuple<typename dof_info<R>::dof_tuple_t,typename dof_info<R>::dof_tuple_t>;
		template<size_t N>
		using dof_t              = std::enable_if_t<(N < count), typename dof_info<R>::template dof_t<N%dof_info<R>::count> >;
		template<size_t N>
		using elem_t             = std::enable_if_t<(N < 2), R>;
		
		template<size_t N>
		static dof_t <N> &dof (value_t &v)    {return (N < count/2) ? dof_get<N>(v.real()) : dof_get<N-count/2>(v.imag());}
		template<size_t N>
		static elem_t<N> &elem(value_t &v)    {return (N==0) ? v.real() : v.imag();}
	};

	template<class Head> // Single-element tuple
	struct dof_info<std::tuple<Head> >
	{
	public:
		// Implementation
		using _iHead = dof_info<Head>;
		template<size_t N> constexpr static bool   _in_head   = (N < _iHead::count);
		template<size_t N> constexpr static size_t _element   = 0;
		template<size_t N> constexpr static size_t _subdegree = N;
	
	public:
		constexpr static size_t
			count                = dof_info<Head>::count,
			elems                = 1;
		constexpr static bool
			primitive            = false,
			primitive_continuous = false,
			primitive_discrete   = false;
		
		using value_t            = std::tuple<Head>;
		using tuple_t            = value_t;
		using dof_tuple_t        = typename dof_info<Head>::dof_tuple_t;
		template<size_t N>
		using dof_t              = std::enable_if_t<_in_head<N>, typename _iHead::template dof_t<N> >;
		template<size_t N>
		using elem_t             = typename std::tuple_element<N, value_t>::type;
		
		template<size_t N>
		static dof_t<N>  &dof (value_t &v)    {return _iHead::dof<N>(std::get<0>(v));}
		template<size_t N>
		static elem_t<N> &elem(value_t &v)    {return std::get<N>(v);}
	};

	template<class Head, class... Tail> // Multi-element tuple
	struct dof_info<std::tuple<Head, Tail...>>
	{
	public:
		// Implementation
		using _iTail = dof_info<std::tuple<Tail...>>;
		using _iHead = dof_info<Head>;
		template<size_t N> constexpr static bool   _in_head = (N < _iHead::count);
		template<size_t N> constexpr static size_t _element = _in_head<N> ? 0 : (1 + _iTail::template _element  <N - _iHead::count>);
		template<size_t N> constexpr static size_t _subdegree = _in_head<N> ? N :    _iTail::template _subdegree<N - _iHead::count>;
	
	public:
		constexpr static size_t
			count                = _iHead::count + _iTail::count,
			elems                = std::tuple_size<std::tuple<Head, Tail...>>::value;
		constexpr static bool
			primitive            = false,
			primitive_continuous = false,
			primitive_discrete   = false;
		
		using value_t            = std::tuple<Head, Tail...>;
		using tuple_t            = value_t;
		using dof_tuple_t        = decltype(std::make_tuple(std::declval<typename _iHead::dof_tuple_t>(), std::declval<typename _iTail::dof_tuple_t>()));
		template<size_t N>
		using dof_t              = std::enable_if_t<(N < count), std::conditional_t<_in_head<N>, typename _iHead::template dof_t<(N)>, typename _iTail::template dof_t<(N - _iHead::count)> > >;
		template<size_t N>
		using elem_t             = typename std::tuple_element<N, value_t>::type;
		
		template<size_t N>
		static dof_t<N>  &dof (value_t &v)    {return dof_get<_subdegree<N>>(std::get<_element<N>>(v));}
		template<size_t N>
		static elem_t<N> &elem(value_t &v)    {return std::get<N>(v);}
	};
	
	// Degrees of freedom for const and volatile types, and references
	#define DOF_Info_CV(QUALIFIERS, SUFFIX, TYPE_MODIFIER) \
		template<class T> \
		struct dof_info<QUALIFIERS T SUFFIX> : public dof_info<T> \
		{ \
			template<size_t N> using  dof_t                               = typename std::TYPE_MODIFIER<typename dof_info<T>::template dof_t <N>>::type; \
			template<size_t N> static dof_t <N> &dof(QUALIFIERS T &v)     {return dof_info<T>::template dof <N>(const_cast<T&>(v));} \
			template<size_t N> using  elem_t                              = typename std::TYPE_MODIFIER<typename dof_info<T>::template elem_t<N>>::type; \
			template<size_t N> static elem_t<N> &elem(QUALIFIERS T &v)    {return dof_info<T>::template elem<N>(const_cast<T&>(v));} \
		}
	DOF_Info_CV(const         , , add_const);
	DOF_Info_CV(      volatile, , add_volatile);
	DOF_Info_CV(const volatile, , add_cv);
	DOF_Info_CV(,              &, add_lvalue_reference);
	DOF_Info_CV(,             &&, add_rvalue_reference);
	#undef DOF_Info_CV
}
