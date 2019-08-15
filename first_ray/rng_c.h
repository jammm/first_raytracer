#ifndef RNG_C
#define RNG_C

// Compile-time RNG taken from https://stackoverflow.com/a/38656025

#include <cstdint>
#include <cstddef>

// Seed
#define RNG_SEED ((__TIME__[7] - '0') * 1  + (__TIME__[6] - '0') * 10  + \
              (__TIME__[4] - '0') * 60   + (__TIME__[3] - '0') * 600 + \
              (__TIME__[1] - '0') * 3600 + (__TIME__[0] - '0') * 36000) + \
              (__LINE__ * 100000)

namespace PRT
{

	typedef uint32_t u32;
	typedef uint64_t u64;
	typedef unsigned char uchar;

	template<u32 S, u32 A = 16807UL, u32 C = 0UL, u32 M = (1UL << 31) - 1>
	struct LinearGenerator {
		static const u32 state = ((u64)S * A + C) % M;
		static const u32 value = state;
		typedef LinearGenerator<state> next;
		struct Split { // Leapfrog
			typedef LinearGenerator< state, A* A, 0, M> Gen1;
			typedef LinearGenerator<next::state, A* A, 0, M> Gen2;
		};
	};

	// Metafunction to get a particular index from generator
	template<u32 S, std::size_t index>
	struct Generate {
		static const uchar value = Generate<LinearGenerator<S>::state, index - 1>::value;
	};

	template<u32 S>
	struct Generate<S, 0> {
		static const uchar value = static_cast<uchar> (LinearGenerator<S>::value);
	};

}


#endif