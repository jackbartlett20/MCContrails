// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#pragma once
#include <array>
#include <algorithm>
#include <numeric>
#include <cmath>

    namespace
calculisto::array
{
/** An array with arithmetic operations.
 * This is a std::array with operators and math functions.
 * (A bit like what std::valarray is to std::vector.)
 */
    template <class T, std::size_t N>
    struct
array_t
{
public:
        using
    backend_type = std::array <T, N>;

        backend_type
    content;

        using
    value_type             = T;
        using
    pointer                = value_type*;
        using
    const_pointer          = const value_type*;
        using
    reference              = value_type&;
        using
    const_reference        = const value_type&;
        using
    iterator               = value_type*;
        using
    const_iterator         = const value_type*;
        using
    size_type              = std::size_t;
        using
    difference_type        = std::ptrdiff_t;
        using
    reverse_iterator       = std::reverse_iterator <iterator>;
        using
    const_reverse_iterator = std::reverse_iterator <const_iterator>;

    array_t () = default;

    /**
     */
    // Converting ctor
        template <
              class... ArgsType
            , class = std::enable_if_t <sizeof...(ArgsType) == N>
        >
        constexpr
    array_t (ArgsType&&... args)
        : content {{ std::forward <ArgsType> (args)... }}
    {}

    // Converting copy ctor
        template <
              class U
            , class = std::enable_if_t <std::is_convertible_v <U, T>>
        >
        explicit (false)
        constexpr
    array_t (array_t <U, N> const& a)
        : content {}
    {
        std::copy (a.begin (), a.end (), begin ());
    }
    /*
        template <class U>
        constexpr array // ditto
    operator= (array <U, N> const& a)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            data_m [i] = a[i];
        }
        //std::copy (a.begin (), a.end (), begin ());
    }
    */
        constexpr void
    fill (const value_type& u)
    {
        content.fill (u);
    }
    // Iterators.
        constexpr decltype (auto)
    begin () noexcept
    {
       return content.begin ();
    }

        constexpr decltype (auto)
    begin () const noexcept
    {
       return content.begin ();
    }

        constexpr decltype (auto)
    end () noexcept
    {
        return content.end ();
    }

        constexpr decltype (auto)
    end () const noexcept
    {
        return content.end ();
    }

        constexpr decltype (auto)
    rbegin () noexcept
    {
        return content.rbegin ();
    }

        constexpr decltype (auto)
    rbegin () const noexcept
    {
        return content.rbegin ();
    }

        constexpr decltype (auto)
    rend () noexcept
    {
        return content.rend ();
    }

        constexpr decltype (auto)
    rend () const noexcept
    {
        return content.rend ();
    }

        constexpr decltype (auto)
    cbegin () const noexcept
    {
        return content.cbegin ();
    }

        constexpr decltype (auto)
    cend () const noexcept
    {
        return content.cend ();
    }

        constexpr decltype (auto)
    crbegin () const noexcept
    {
        return content.crbegin ();
    }

        constexpr decltype (auto)
    crend () const noexcept
    {
        return content.crend ();
    }

    // Capacity.
        constexpr size_type
    size () const noexcept
    {
        return content.size ();
    }

        constexpr size_type
    max_size () const noexcept
    {
        return content.max_size ();
    }

        constexpr bool
    empty () const noexcept
    {
        return content.empty ();
    }

    // Element access.
        constexpr decltype (auto)
    operator[] (size_type n) noexcept
    {
        return content[n];
    }

        constexpr decltype (auto)
    operator[] (size_type n) const noexcept
    {
        return content[n];
    }

        constexpr decltype (auto)
    at (size_type n)
    {
        return content.at (n);
    }

        constexpr decltype (auto)
    at (size_type n) const
    {
        return content.at (n);
    }

        constexpr decltype (auto)
    front () noexcept
    {
        return content.front ();
    }

        constexpr decltype (auto)
    front () const noexcept
    {
        return content.front ();
    }

        constexpr decltype (auto)
    back () noexcept
    {
        return content.back ();
    }

        constexpr decltype (auto)
    back () const noexcept
    {
        return content.back ();
    }
        constexpr decltype (auto)
    data () noexcept
    {
        return content.data ();
    }
        constexpr decltype (auto)
    data () const noexcept
    {
        return content.data ();
    }

        /*
    // ---
    // unary operators:
        constexpr array
    operator+ () const
    {
            auto
        r = *this;
        std::transform (r.begin (), r.end (), r.begin (), [&](auto e){ return +e; });
        return r;
    }
        constexpr array
    operator- () const
    {
            auto
        r = *this;
        std::transform (r.begin (), r.end (), r.begin (), [&](auto e){ return -e; });
        return r;
    }
    */


#   define COMPUTED_ASSIGNEMENT(op)               \
        constexpr auto                            \
    operator op##= (const value_type& x)          \
    {                                             \
        std::transform (                          \
                begin (),                         \
                end (),                           \
                begin (),                         \
                [=](auto& e){ return e op##= x; } \
        );                                        \
        return *this;                             \
    }
    COMPUTED_ASSIGNEMENT(*)
    COMPUTED_ASSIGNEMENT(/)
    COMPUTED_ASSIGNEMENT(%)
    COMPUTED_ASSIGNEMENT(+)
    COMPUTED_ASSIGNEMENT(-)
    COMPUTED_ASSIGNEMENT(^)
    COMPUTED_ASSIGNEMENT(&)
    COMPUTED_ASSIGNEMENT(|)
    COMPUTED_ASSIGNEMENT(<<)
    COMPUTED_ASSIGNEMENT(>>)
#   undef COMPUTED_ASSIGNEMENT
        /*
#   define COMPUTED_ASSIGNEMENT(op)                                                      \
        template <                                                                       \
              class U                                                                    \
            , class = decltype (std::declval <value_type> () op##= std::declval <U> ()) \
        >                                                                                \
        constexpr array&                                                                 \
    operator op##= (U const& x)                                                          \
    {                                                                                    \
        for (auto& e: data_m)                                                            \
        {                                                                                \
            e op##= x;                                                                   \
        }                                                                                \
        return *this;                                                                    \
    }
    COMPUTED_ASSIGNEMENT(*)
    COMPUTED_ASSIGNEMENT(/)
    COMPUTED_ASSIGNEMENT(%)
    COMPUTED_ASSIGNEMENT(+)
    COMPUTED_ASSIGNEMENT(-)
    COMPUTED_ASSIGNEMENT(^)
    COMPUTED_ASSIGNEMENT(&)
    COMPUTED_ASSIGNEMENT(|)
    COMPUTED_ASSIGNEMENT(<<)
    COMPUTED_ASSIGNEMENT(>>)
#   undef COMPUTED_ASSIGNEMENT
*/
// Don't know how to constraint that...
// , class = decltype (std::declval <value_type> () op##= std::declval <U> ())
#   define COMPUTED_ASSIGNEMENT(op)                                                      \
        template <                                                                       \
              class U                                                                    \
        >                                                                                \
        constexpr decltype (auto)                                                        \
    operator op##= (array_t <U, N> const& v)                                             \
    {                                                                                    \
        std::transform (                                                                 \
            begin   (),                                                                  \
            end     (),                                                                  \
            v.begin (),                                                                  \
            begin   (),                                                                  \
            [=](auto a, auto b){ return a op b; }                                        \
        );                                                                               \
        return *this;                                                                    \
    }
    COMPUTED_ASSIGNEMENT(*)
    COMPUTED_ASSIGNEMENT(/)
    COMPUTED_ASSIGNEMENT(%)
    COMPUTED_ASSIGNEMENT(+)
    COMPUTED_ASSIGNEMENT(-)
    COMPUTED_ASSIGNEMENT(^)
    COMPUTED_ASSIGNEMENT(&)
    COMPUTED_ASSIGNEMENT(|)
    COMPUTED_ASSIGNEMENT(<<)
    COMPUTED_ASSIGNEMENT(>>)
#   undef COMPUTED_ASSIGNEMENT

        auto
    operator <=> (array_t <T, N> const&) const = default;
}; // class array

// TODO swap

    template <class T, class... U>
array_t (T, U...) -> array_t <T, 1 + sizeof... (U)>;
    template <
          class T
        , std::size_t N
        , class V = decltype (+std::declval <T> ())
    >
    constexpr array_t <V, N>
operator + (array_t <T, N> const& a)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          a.begin ()
        , a.end   ()
        , r.begin ()
        , [](auto x){ return +x; }
    );
    return r;
}
    template <
          class T
        , std::size_t N
        , class V = decltype (-std::declval <T> ())
    >
    constexpr array_t <V, N>
operator - (array_t <T, N> const& a)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          a.begin ()
        , a.end   ()
        , r.begin ()
        , [](auto x){ return -x; }
    );
    return r;
}

    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () + std::declval <U> ())
    >
    constexpr array_t <V, N>
operator + (array_t <T, N> const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          a.begin ()
        , a.end   ()
        , b.begin ()
        , r.begin ()
        , [](auto x, auto y){ return x + y; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () - std::declval <U> ())
    >
    constexpr array_t <V, N>
operator - (array_t <T, N> const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          a.begin ()
        , a.end   ()
        , b.begin ()
        , r.begin ()
        , [](auto x, auto y){ return x - y; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () * std::declval <U> ())
    >
    constexpr array_t <V, N>
operator * (array_t <T, N> const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          a.begin ()
        , a.end   ()
        , b.begin ()
        , r.begin ()
        , [](auto x, auto y){ return x * y; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () / std::declval <U> ())
    >
    constexpr array_t <V, N>
operator / (array_t <T, N> const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          a.begin ()
        , a.end   ()
        , b.begin ()
        , r.begin ()
        , [](auto x, auto y){ return x / y; }
    );
    return r;
}

    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () + std::declval <U> ())
    >
    constexpr array_t <V, N>
operator + (array_t <T, N> const& a, U const& b)
{
        auto
    r = array_t <V, N> {};
    // TODO capture by value or by ref depending on U (using if constexpr)?
    // E.g. is U is trivially copyable?
    std::transform (
          a.begin ()
        , a.end   ()
        , r.begin ()
        , [&](auto x){ return x + b; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () + std::declval <U> ())
    >
    constexpr array_t <V, N>
operator + (T const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          b.begin ()
        , b.end   ()
        , r.begin ()
        , [&](auto x){ return a + x; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () - std::declval <U> ())
    >
    constexpr array_t <V, N>
operator - (array_t <T, N> const& a, U const& b)
{
        auto
    r = array_t <V, N> {};
    // TODO capture by value or by ref depending on U (using if constexpr)?
    // E.g. is U is trivially copyable?
    std::transform (
          a.begin ()
        , a.end   ()
        , r.begin ()
        , [&](auto x){ return x - b; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () - std::declval <U> ())
    >
    constexpr array_t <V, N>
operator - (T const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          b.begin ()
        , b.end   ()
        , r.begin ()
        , [&](auto x){ return a - x; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () * std::declval <U> ())
    >
    constexpr array_t <V, N>
operator * (array_t <T, N> const& a, U const& b)
{
        auto
    r = array_t <V, N> {};
    // TODO capture by value or by ref depending on U (using if constexpr)?
    // E.g. is U is trivially copyable?
    std::transform (
          a.begin ()
        , a.end   ()
        , r.begin ()
        , [&](auto x){ return x * b; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () * std::declval <U> ())
    >
    constexpr array_t <V, N>
operator * (T const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          b.begin ()
        , b.end   ()
        , r.begin ()
        , [&](auto x){ return a * x; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () / std::declval <U> ())
    >
    constexpr array_t <V, N>
operator / (array_t <T, N> const& a, U const& b)
{
        auto
    r = array_t <V, N> {};
    // TODO capture by value or by ref depending on U (using if constexpr)?
    // E.g. is U is trivially copyable?
    std::transform (
          a.begin ()
        , a.end   ()
        , r.begin ()
        , [&](auto x){ return x / b; }
    );
    return r;
}
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (std::declval <T> () / std::declval <U> ())
    >
    constexpr array_t <V, N>
operator / (T const& a, array_t <U, N> const& b)
{
        auto
    r = array_t <V, N> {};
    std::transform (
          b.begin ()
        , b.end   ()
        , r.begin ()
        , [&](auto x){ return a / x; }
    );
    return r;
}

#define UNARY_FUNCTION(f)                              \
    template <                                         \
          class T                                      \
        , std::size_t N                                \
        , class U = decltype (f (std::declval <T> ())) \
    >                                                  \
    constexpr auto                                     \
f (array_t <T, N> const& p)                            \
{                                                      \
        auto                                           \
    r = array_t <U, N> {};                             \
    using std::f;                                      \
    std::transform (                                   \
            p.begin (),                                \
            p.end   (),                                \
            r.begin (),                                \
            [](auto e){ return f (e); }                \
    );                                                 \
    return r;                                          \
}
// See below.
// UNARY_FUNCTION(abs)
UNARY_FUNCTION(cos)
UNARY_FUNCTION(acos)
UNARY_FUNCTION(cosh)
UNARY_FUNCTION(sin)
UNARY_FUNCTION(asin)
UNARY_FUNCTION(sinh)
UNARY_FUNCTION(tan)
UNARY_FUNCTION(tanh)
UNARY_FUNCTION(atan)
// UNARY_FUNCTION(atan2)
UNARY_FUNCTION(exp)
UNARY_FUNCTION(log)
UNARY_FUNCTION(log10)
UNARY_FUNCTION(sqrt)
#undef UNARY_FUNCTION

/* Can't use the above macro because
 * both std::abs and ::abs exist and this make
 * the return type deduction fail.
 *
 */
    template <
          class T
        , std::size_t N
    > constexpr auto
abs (array_t <T, N> p)
{
    using std::abs;
    std::transform (
            p.begin (),
            p.end   (),
            p.begin (),
            [=](auto e){ return abs (e); }
    );
    return p;
}
    using std::pow; // <- TODO: Is there a way to remove that?
    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (pow (std::declval <T> (), std::declval <U> ()))
    >
    constexpr array_t <V, N>
pow (array_t <T, N> const& p, U exp)
{
        auto
    r = array_t <V, N> {};
    std::transform (
            p.begin (),
            p.end   (),
            r.begin (),
            [=](auto e){ return pow (e, exp); }
    );
    return r;
}

    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (pow (std::declval <T> (), std::declval <U> ()))
    >
    constexpr array_t <V, N>
pow (T t, array_t <U, N> const& exp)
{
        auto
    r = array_t <V, N> {};
    std::transform (
            exp.begin (),
            exp.end   (),
            r.begin (),
            [=](auto ee){ return pow (t, ee); }
    );
    return r;
}

    template <
          class T
        , class U
        , std::size_t N
        , class V = decltype (pow (std::declval <T> (), std::declval <U> ()))
    >
    constexpr array_t <V, N>
pow (array_t <T, N> const& p, array_t <U, N> const& exp)
{
        auto
    r = array_t <V, N> {};
    std::transform (
            p.begin (),
            p.end   (),
            exp.begin (),
            r.begin (),
            [=](auto e, auto ee){ return pow (e, ee); }
    );
    return r;
}
    template <class T, std::size_t N>
    constexpr auto
hypot (array_t <T, N> const& p)
{
        using std::sqrt;
    p = p * p;
    return sqrt (std::accumulate (
        std::begin (p),
        std::end   (p),
        static_cast <T> (0)
    ));
}
    template <class T>
    constexpr auto
hypot (array_t <T, 2> const& p)
{
        using std::hypot;
    return hypot (p[0], p[1]);
}
    template <class T>
    constexpr auto
hypot (array_t <T, 3> const& p)
{
        using std::hypot;
    return hypot (p[0], p[1], p[2]);
}

// Concatenation.
    namespace
detail
{
        template <
              class T
            , std::size_t    N1
            , std::size_t... I1
            , std::size_t    N2
            , std::size_t... I2
        >
        constexpr auto
    concatenate (
              array_t <T, N1> const& a1
            , array_t <T, N2> const& a2
            , std::index_sequence <I1...>
            , std::index_sequence <I2...>
    ){
        return array_t <T, N1 + N2> { a1[I1]..., a2[I2]... };
    }
} // namespace detail
    template <
          class T
        , std::size_t N1
        , std::size_t N2
    >
    constexpr array_t <T, N1 + N2>
concatenate (
          array_t <T, N1> const& a1
        , array_t <T, N2> const& a2
){
    return detail::concatenate (a1, a2, std::make_index_sequence <N1> (), std::make_index_sequence <N2> ());
}

    template <class T, std::size_t N1, std::size_t N2>
    array_t <T, N1 + N2>
operator , (array_t <T, N1> const& a, array_t <T, N2> const& b)
{
    return concatenate (a, b);
}
    template <class... Arrays>
    auto
concatenate (Arrays&&... arrays)
{
    return ( ... ,  arrays);
}

    template <class T, auto N> 
    constexpr auto
sum (array_t <T, N> const& array)
{
    return std::accumulate (
          array.begin ()
        , array.end ()
        , static_cast <T> (0)
    );
}


} // namespace calculisto::array
