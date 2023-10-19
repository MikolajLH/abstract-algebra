#pragma once


namespace prime
{
	//https://en.wikipedia.org/wiki/Integer_square_root
	constexpr size_t int_sqrt_bin_search(size_t n)
	{
		size_t L = 0;
		size_t R = n + 1;
		while (L != R - 1)
		{
			const size_t M = (L + R) / 2;
			if (M * M <= n)
				L = M;
			else
				R = M;
		}
		return L;
	}

	//https://en.wikipedia.org/wiki/Integer_square_root
	constexpr size_t int_sqrt_heron(size_t n)
	{
		if (n <= 1)return n;
		size_t x0 = n / 2;
		size_t x1 = (x0 + n / x0) / 2;
		while (x1 < x0)
		{
			x0 = x1;
			x1 = (x0 + n / x0) / 2;
		}
		return x0;
	}

	namespace impl
	{
		template<size_t P, size_t Q>
		struct IsPrimeImpl { inline static constexpr bool value = P % Q != 0 and IsPrimeImpl<P, Q - 1>::value; };

		template<size_t P>
		struct IsPrimeImpl<P, 2> { inline static constexpr bool value = P % 2 != 0; };

		template<size_t P>
		struct IsPrimeStruct { inline static constexpr bool value = IsPrimeImpl<P, int_sqrt_heron(P)>::value; };

		template<> struct IsPrimeStruct<0> { inline static constexpr bool value = false; };
		template<> struct IsPrimeStruct<1> { inline static constexpr bool value = false; };

	}

	template<size_t P>
	inline constexpr bool IsPrime = impl::IsPrimeStruct<P>::value;


	constexpr bool is_prime(size_t n)
	{
		return false;
	}

#define PRIME(P) template<> struct prime::impl::IsPrimeStruct<P> { inline static constexpr bool value = true; }

	PRIME(2);    PRIME(3);    PRIME(5);    PRIME(7);    PRIME(11);   PRIME(13);   PRIME(17);   PRIME(19);   PRIME(23);   PRIME(29);
	PRIME(31);   PRIME(37);   PRIME(41);   PRIME(43);   PRIME(47);   PRIME(53);   PRIME(59);   PRIME(61);   PRIME(67);   PRIME(71);
	PRIME(73);   PRIME(79);   PRIME(83);   PRIME(89);   PRIME(97);   PRIME(101);  PRIME(103);  PRIME(107);  PRIME(109);  PRIME(113);
	PRIME(127);  PRIME(131);  PRIME(137);  PRIME(139);  PRIME(149);  PRIME(151);  PRIME(157);  PRIME(163);  PRIME(167);  PRIME(173);
	PRIME(179);  PRIME(181);  PRIME(191);  PRIME(193);  PRIME(197);  PRIME(199);  PRIME(211);  PRIME(223);  PRIME(227);  PRIME(229);
	PRIME(233);  PRIME(239);  PRIME(241);  PRIME(251);  PRIME(257);  PRIME(263);  PRIME(269);  PRIME(271);  PRIME(277);  PRIME(281);
	PRIME(283);  PRIME(293);  PRIME(307);  PRIME(311);  PRIME(313);  PRIME(317);  PRIME(331);  PRIME(337);  PRIME(347);  PRIME(349);
	PRIME(353);  PRIME(359);  PRIME(367);  PRIME(373);  PRIME(379);  PRIME(383);  PRIME(389);  PRIME(397);  PRIME(401);  PRIME(409);
	PRIME(419);  PRIME(421);  PRIME(431);  PRIME(433);  PRIME(439);  PRIME(443);  PRIME(449);  PRIME(457);  PRIME(461);  PRIME(463);
	PRIME(467);  PRIME(479);  PRIME(487);  PRIME(491);  PRIME(499);  PRIME(503);  PRIME(509);  PRIME(521);  PRIME(523);  PRIME(541);
	PRIME(547);  PRIME(557);  PRIME(563);  PRIME(569);  PRIME(571);  PRIME(577);  PRIME(587);  PRIME(593);  PRIME(599);  PRIME(601);
	PRIME(607);  PRIME(613);  PRIME(617);  PRIME(619);  PRIME(631);  PRIME(641);  PRIME(643);  PRIME(647);  PRIME(653);  PRIME(659);
	PRIME(661);  PRIME(673);  PRIME(677);  PRIME(683);  PRIME(691);  PRIME(701);  PRIME(709);  PRIME(719);  PRIME(727);  PRIME(733);
	PRIME(739);  PRIME(743);  PRIME(751);  PRIME(757);  PRIME(761);  PRIME(769);  PRIME(773);  PRIME(787);  PRIME(797);  PRIME(809);
	PRIME(811);  PRIME(821);  PRIME(823);  PRIME(827);  PRIME(829);  PRIME(839);  PRIME(853);  PRIME(857);  PRIME(859);  PRIME(863);
	PRIME(877);  PRIME(881);  PRIME(883);  PRIME(887);  PRIME(907);  PRIME(911);  PRIME(919);  PRIME(929);  PRIME(937);  PRIME(941);
	PRIME(947);  PRIME(953);  PRIME(967);  PRIME(971);  PRIME(977);  PRIME(983);  PRIME(991);  PRIME(997);  PRIME(1009); PRIME(1013);
	PRIME(1019); PRIME(1021); PRIME(1031); PRIME(1033); PRIME(1039); PRIME(1049); PRIME(1051); PRIME(1061); PRIME(1063); PRIME(1069);
	PRIME(1087); PRIME(1091); PRIME(1093); PRIME(1097); PRIME(1103); PRIME(1109); PRIME(1117); PRIME(1123); PRIME(1129); PRIME(1151);
	PRIME(1153); PRIME(1163); PRIME(1171); PRIME(1181); PRIME(1187); PRIME(1193); PRIME(1201); PRIME(1213); PRIME(1217); PRIME(1223);
}