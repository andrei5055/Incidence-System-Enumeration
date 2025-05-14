#pragma once
#include "Global.h"

typedef unsigned int uint;
#define E  0x80000000

class Combination {
public:
	Combination(tchar a, tchar b, tchar c) : m_code((E>>a) + (E>>b) + (E>>c)) {}
	inline auto code() const		{ return m_code; }
private:
	const uint m_code;
	Combination* m_pPrev;
	Combination* m_pNext;
};

class CombinationSet {
public:
	CombinationSet(int n, int k);
};
