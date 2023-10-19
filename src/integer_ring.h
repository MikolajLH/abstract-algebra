#pragma once
#include <compare>
#include "primes.h"
#include <print>
#include <format>
#include <concepts>
#include <cassert>
#include <numeric>
#include <set>
#include <math.h>
#include <limits>

namespace algebra
{
	template<size_t N>
	class Zn
	{
		public:
			Zn() = delete;
		
			constexpr Zn(std::unsigned_integral auto val) noexcept 
				: m_val{ size_t(val) % N } 
			{}

			constexpr Zn(std::signed_integral auto val) noexcept 
				: m_val{ size_t(val < 0 ? -val : val) % N } 
			{}

			inline constexpr auto operator <=>(const Zn<N>&) const noexcept = default;

			inline constexpr auto operator +(this const auto& self, const Zn<N>& other) noexcept
			{
				if constexpr (N - 1 > std::numeric_limits<size_t>::max() / 2)
				{
					const size_t s = self.integer();
					const size_t o = other.integer();

					const size_t minus_s = N - self.integer();
					const size_t minus_o = N - other.integer();

					return Zn<N>(s - (s > minus_o) * minus_o + (not (s > minus_o)) * o);

				}
				else return Zn<N>(self.m_val + other.m_val);
			}

			inline constexpr auto operator *(this const auto& self, const Zn<N>& other) noexcept
			{
				static_assert(N - 1 <= prime::int_sqrt_heron(std::numeric_limits<size_t>::max()), "Possible overflow");

				return Zn<N>(self.m_val * other.m_val);
			}

			inline constexpr void operator++(this auto& self) noexcept { self = self + 1; }

			inline constexpr void operator--(this auto& self) noexcept { self = self - 1; }

			inline constexpr auto operator -(this const auto& self) noexcept
			{
				return Zn<N>(N - self.m_val);
			}

			inline constexpr auto operator -(this const auto& self, const Zn<N>& other) noexcept
			{
				return self + (-other);
			}
			
			constexpr auto inverse(this const auto& self)noexcept requires prime::IsPrime<N>
			{
				assert(self.integer() != 0);
				size_t r0 = self.integer();
				size_t r1 = self.modulus();
				Zn<N> s0 = 1;
				Zn<N> s1 = 0;

				while (r1 != 0)
				{
					const size_t qi = r0 / r1;
					r0 = std::exchange(r1, r0 % r1);
					s0 = std::exchange(s1, s0 - s1 * qi);
				}

				return s0;
			}
			
			constexpr auto operator /(this const auto& self, const Zn<N>& other) requires prime::IsPrime<N>
			{
				static_assert(other.integer() != 0, "Cannot divide by zero!");
				assert(other.integer() == 0);
				return self * other.inverse();
			}

			static constexpr size_t modulus() noexcept { return N; }

			inline constexpr size_t integer() const noexcept { return m_val; }

			static constexpr Zn<N> primitive_root() noexcept
			{
				if constexpr (prime::IsPrime<N>)
				{
					
					for (size_t m = 0; m < N; ++m)
					{
						bool remainders[N] = {};

						Zn<N> k = m;
						for (size_t i = 0; i < N; ++i, k = k * m)
						{
							remainders[k.integer()] = true;
						}

						if (not remainders[0] and std::accumulate(remainders + 1, remainders + N, true, [](bool a, bool b) {return a and b; }))
							return m;
					}

					//there always exists primitive root if N is a prime number
					std::unreachable();
				}
				else
				{
					/*TODO for other cases*/
					static_assert(prime::IsPrime<N>, "N has to be prime");
				}
			}

			static constexpr std::set<Zn<N>> primitive_roots() noexcept
			{
				if constexpr (prime::IsPrime<N>)
				{
					std::set<Zn<N>> roots;

					for (size_t m = 0; m < N; ++m)
					{
						bool remainders[N] = {};

						Zn<N> k = m;
						for (size_t i = 0; i < N; ++i, k = k * m)
						{
							remainders[k.integer()] = true;
						}

						if (not remainders[0] and std::accumulate(remainders + 1, remainders + N, true, [](bool a, bool b) {return a and b; }))
							roots.emplace(m);
					}

					return roots;
				}
				else
				{
					/*TODO for other cases*/
					static_assert(prime::IsPrime<N>, "N has to be prime");
				}
			}

		private:
			size_t m_val;
	};

#define DECLARE_Zn_literal(N) constexpr auto operator""_Z##N(size_t val){return algebra::Zn<N>(val);}
#define DECLARE_Zn_conv_fun(N) constexpr auto Z##N(std::integral auto val){return algebra::Zn<N>(val);}
}

template<size_t N>
struct std::formatter<algebra::Zn<N>> {
	constexpr auto parse(std::format_parse_context& ctx) {
		return ctx.begin();
	}

	auto format(const algebra::Zn<N>& obj, std::format_context& ctx) {
		return std::format_to(ctx.out(), "{}_Z{}", obj.integer(), obj.modulus());
	}
};


//https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
//returns (gcd(a,b), s, t) such that a * s + b * t = gcd(a,b)
template<class Int>
[[nodiscard("pure function")]] constexpr auto extended_euclidean_algorithm(Int a, Int b)
{
	static_assert(not std::unsigned_integral<Int>, "since this function uses substraction, using unsigned integral types can couse underflow");

	Int r0 = std::max(a, b);
	Int r1 = std::min(a, b);

	Int s0 = 1;
	Int s1 = 0;

	Int t0 = 0;
	Int t1 = 1;

	while (r1 != 0)
	{
		const Int qi = r0 / r1;
		r0 = std::exchange(r1, r0 - qi * r1);
		s0 = std::exchange(s1, s0 - qi * s1);
		t0 = std::exchange(t1, t0 - qi * t1);
	}

	return std::make_tuple(r0, s0, t0);
}