#ifndef NAGISS_LIBRARY_HPP
#define NAGISS_LIBRARY_HPP

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include<iostream>
#include<iomanip>
#include<vector>
#include<set>
#include<map>
#include<unordered_set>
#include<unordered_map>
#include<algorithm>
#include<numeric>
#include<limits>
#include<bitset>
#include<functional>
#include<type_traits>
#include<queue>
#include<stack>
#include<array>
#include<random>
#include<utility>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<ctime>
#include<string>
#include<sstream>
#include<chrono>
#include<climits>

#ifdef _MSC_VER
#include<intrin.h>
#endif
#ifdef __GNUC__
#include<x86intrin.h>
#endif

#ifdef __GNUC__
//#pragma GCC target("avx2")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native")
#pragma GCC optimize("O3")
#pragma GCC optimize("Ofast")
//#pragma GCC optimize("unroll-loops")

#pragma clang attribute push (__attribute__((target("arch=skylake"))),apply_to=function)
/* 最後に↓を貼る
#ifdef __GNUC__
#pragma clang attribute pop
#endif
*/
#endif


// ========================== macroes ==========================

#define rep(i,n) for(ll i=0; (i)<(n); (i)++)
#define rep1(i,n) for(ll i=1; (i)<=(n); (i)++)
#define rep3(i,l,r) for(auto i=(l); (i)<(r); (i)++)

//#define NDEBUG

#ifndef NDEBUG
#define ASSERT(expr, ...) \
		do { \
			if(!(expr)){ \
				printf("%s(%d): Assertion failed.\n", __FILE__, __LINE__); \
				printf(__VA_ARGS__); \
				abort(); \
			} \
		} while (false)
#else
#define ASSERT(...)
#endif

#define ASSERT_RANGE(value, left, right) \
	ASSERT((left <= value) && (value < right), \
		"`%s` (%d) is out of range [%d, %d)", #value, value, left, right)

#define CHECK(var) do{ std::cout << #var << '=' << var << endl; } while (false)

// ========================== utils ==========================

using namespace std;
using ll = long long;
constexpr double PI = 3.1415926535897932;

template<class T, class S> inline bool chmin(T& m, S q) {
	if (m > q) { m = q; return true; }
	else return false;
}

template<class T, class S> inline bool chmax(T& m, const S q) {
	if (m < q) { m = q; return true; }
	else return false;
}

// クリッピング  // clamp (C++17) と等価
template<class T> inline T clipped(const T& v, const T& low, const T& high) {
	return min(max(v, low), high);
}

// 2 次元ベクトル
template<typename T> struct Vec2 {
	/*
	y 軸正は下方向
	x 軸正は右方向
	回転は時計回りが正（y 軸正を上と考えると反時計回りになる）
	*/
	T y, x;
	constexpr inline Vec2() = default;
	constexpr inline Vec2(const T& arg_y, const T& arg_x) : y(arg_y), x(arg_x) {}
	inline Vec2(const Vec2&) = default;  // コピー
	inline Vec2(Vec2&&) = default;  // ムーブ
	inline Vec2& operator=(const Vec2&) = default;  // 代入
	inline Vec2& operator=(Vec2&&) = default;  // ムーブ代入
	template<typename S> constexpr inline Vec2(const Vec2<S>& v) : y((T)v.y), x((T)v.x) {}
	inline Vec2 operator+(const Vec2& rhs) const {
		return Vec2(y + rhs.y, x + rhs.x);
	}
	inline Vec2 operator+(const T& rhs) const {
		return Vec2(y + rhs, x + rhs);
	}
	inline Vec2 operator-(const Vec2& rhs) const {
		return Vec2(y - rhs.y, x - rhs.x);
	}
	template<typename S> inline Vec2 operator*(const S& rhs) const {
		return Vec2(y * rhs, x * rhs);
	}
	inline Vec2 operator*(const Vec2& rhs) const {  // x + yj とみなす
		return Vec2(x * rhs.y + y * rhs.x, x * rhs.x - y * rhs.y);
	}
	template<typename S> inline Vec2 operator/(const S& rhs) const {
		ASSERT(rhs != 0.0, "Zero division!");
		return Vec2(y / rhs, x / rhs);
	}
	inline Vec2 operator/(const Vec2& rhs) const {  // x + yj とみなす
		return (*this) * rhs.inv();
	}
	inline Vec2& operator+=(const Vec2& rhs) {
		y += rhs.y;
		x += rhs.x;
		return *this;
	}
	inline Vec2& operator-=(const Vec2& rhs) {
		y -= rhs.y;
		x -= rhs.x;
		return *this;
	}
	template<typename S> inline Vec2& operator*=(const S& rhs) const {
		y *= rhs;
		x *= rhs;
		return *this;
	}
	inline Vec2& operator*=(const Vec2& rhs) {
		*this = (*this) * rhs;
		return *this;
	}
	inline Vec2& operator/=(const Vec2& rhs) {
		*this = (*this) / rhs;
		return *this;
	}
	inline bool operator!=(const Vec2& rhs) const {
		return x != rhs.x || y != rhs.y;
	}
	inline bool operator==(const Vec2& rhs) const {
		return x == rhs.x && y == rhs.y;
	}
	inline void rotate(const double& rad) {
		*this = rotated(rad);
	}
	inline Vec2<double> rotated(const double& rad) const {
		return (*this) * rotation(rad);
	}
	static inline Vec2<double> rotation(const double& rad) {
		return Vec2(sin(rad), cos(rad));
	}
	static inline Vec2<double> rotation_deg(const double& deg) {
		return rotation(PI * deg / 180.0);
	}
	inline Vec2<double> rounded() const {
		return Vec2<double>(round(y), round(x));
	}
	inline Vec2<double> inv() const {  // x + yj とみなす
		const double norm_sq = l2_norm_square();
		ASSERT(norm_sq != 0.0, "Zero division!");
		return Vec2(-y / norm_sq, x / norm_sq);
	}
	inline double l2_norm() const {
		return sqrt(x * x + y * y);
	}
	inline double l2_norm_square() const {
		return x * x + y * y;
	}
	inline T l1_norm() const {
		return std::abs(x) + std::abs(y);
	}
	inline double abs() const {
		return l2_norm();
	}
	inline double phase() const {  // [-PI, PI) のはず
		return atan2(y, x);
	}
	inline double phase_deg() const {  // [-180, 180) のはず
		return phase() / PI * 180.0;
	}
};
template<typename T, typename S> inline Vec2<T> operator*(const S& lhs, const Vec2<T>& rhs) {
	return rhs * lhs;
}
template<typename T> ostream& operator<<(ostream& os, const Vec2<T>& vec) {
	os << vec.y << ' ' << vec.x;
	return os;
}

// 乱数
struct Random {
	using ull = unsigned long long;
	ull seed;
	inline Random(ull aSeed) : seed(aSeed) {
		ASSERT(seed != 0ull, "Seed should not be 0.");
	}
	const inline ull& next() {
		seed ^= seed << 9;
		seed ^= seed >> 7;
		return seed;
	}
	// (0.0, 1.0)
	inline double random() {
		return (double)next() / (double)ULLONG_MAX;
	}
	// [0, right)
	inline int randint(const int right) {
		return next() % (ull)right;
	}
	// [left, right)
	inline int randint(const int left, const int right) {
		return next() % (ull)(right - left) + left;
	}
};

// 2 次元配列
template<class T, int height, int width> struct Board {
	array<T, height * width> data;
	template<class Int> constexpr inline auto& operator[](const Vec2<Int>& p) {
		return data[width * p.y + p.x];
	}
	template<class Int> constexpr inline const auto& operator[](const Vec2<Int>& p) const {
		return data[width * p.y + p.x];
	}
	template<class Int> constexpr inline auto& operator[](const initializer_list<Int>& p) {
		return data[width * *p.begin() + *(p.begin() + 1)];
	}
	template<class Int> constexpr inline const auto& operator[](const initializer_list<Int>& p) const {
		return data[width * *p.begin() + *(p.begin() + 1)];
	}
	constexpr inline void Fill(const T& fill_value) {
		fill(data.begin(), data.end(), fill_value);
	}
	void Print() const {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				cout << data[width * y + x] << " \n"[x == width - 1];
			}
		}
	}
};

// キュー
template<class T, int max_size> struct Queue {
	array<T, max_size> data;
	int left, right;
	inline Queue() : data(), left(0), right(0) {}
	inline Queue(initializer_list<T> init) :
		data(init.begin(), init.end()), left(0), right(init.size()) {}

	inline bool empty() const {
		return left == right;
	}
	inline void push(const T& value) {
		data[right] = value;
		right++;
	}
	inline void pop() {
		left++;
	}
	const inline T& front() const {
		return data[left];
	}
	template <class... Args> inline void emplace(const Args&... args) {
		data[right] = T(args...);
		right++;
	}
	inline void clear() {
		left = 0;
		right = 0;
	}
	inline int size() const {
		return right - left;
	}
};

// スタック  // コンストラクタ呼ぶタイミングとかが考えられてなくて良くない
template<class T, int max_size> struct Stack {
	array<T, max_size> data;
	int right;
	inline Stack() : data(), right(0) {}
	inline Stack(const int n) : data(), right(0) { resize(n); }
	inline Stack(const int n, const T& val) : data(), right(0) { resize(n, val); }
	inline Stack(const initializer_list<T>& init) : data(), right(init.size()) {
		memcpy(&data[0], init.begin(), sizeof(T) * init.size());
	}  // これ大丈夫か？
	inline Stack(const Stack& rhs) : data(), right(rhs.right) {  // コピー
		for (int i = 0; i < right; i++) {
			data[i] = rhs.data[i];
		}
	}
	template<class S> inline Stack(const Stack<S, max_size>& rhs) : data(), right(rhs.right) {
		for (int i = 0; i < right; i++) {
			data[i] = rhs.data[i];
		}
	}
	Stack& operator=(const Stack& rhs) {
		right = rhs.right;
		for (int i = 0; i < right; i++) {
			data[i] = rhs.data[i];
		}
		return *this;
	}
	Stack& operator=(const vector<T>& rhs) {
		right = (int)rhs.size();
		ASSERT(right <= max_size, "too big vector");
		for (int i = 0; i < right; i++) {
			data[i] = rhs[i];
		}
		return *this;
	}
	Stack& operator=(Stack&&) = default;
	inline bool empty() const {
		return 0 == right;
	}
	inline void push(const T& value) {
		ASSERT_RANGE(right, 0, max_size);
		data[right] = value;
		right++;
	}
	inline T pop() {
		right--;
		ASSERT_RANGE(right, 0, max_size);
		return data[right];
	}
	const inline T& top() const {
		return data[right - 1];
	}
	template <class... Args> inline void emplace(Args&&... args) {
		ASSERT_RANGE(right, 0, max_size);
		new(&data[right])T(forward(args...));
		right++;
	}
	inline void clear() {
		right = 0;
	}
	inline void insert(const int& idx, const T& value) {
		ASSERT_RANGE(idx, 0, right + 1);
		ASSERT_RANGE(right, 0, max_size);
		int i = right;
		right++;
		while (i != idx) {
			data[i] = data[i - 1];
			i--;
		}
		data[idx] = value;
	}
	inline void del(const int& idx) {
		ASSERT_RANGE(idx, 0, right);
		right--;
		for (int i = idx; i < right; i++) {
			data[i] = data[i + 1];
		}
	}
	inline int index(const T& value) const {
		for (int i = 0; i < right; i++) {
			if (value == data[i]) return i;
		}
		return -1;
	}
	inline void remove(const T& value) {
		int idx = index(value);
		ASSERT(idx != -1, "not contain the value.");
		del(idx);
	}
	inline void resize(const int& sz) {
		ASSERT_RANGE(sz, 0, max_size + 1);
		for (; right < sz; right++) {
			data[right].~T();
			new(&data[right]) T();
		}
		right = sz;
	}
	inline void resize(const int& sz, const T& fill_value) {
		ASSERT_RANGE(sz, 0, max_size + 1);
		for (; right < sz; right++) {
			data[right].~T();
			new(&data[right]) T(fill_value);
		}
		right = sz;
	}
	inline int size() const {
		return right;
	}
	inline T& operator[](const int n) {
		ASSERT_RANGE(n, 0, right);
		return data[n];
	}
	inline const T& operator[](const int n) const {
		ASSERT_RANGE(n, 0, right);
		return data[n];
	}
	inline T* begin() {
		return (T*)data.data();
	}
	inline const T* begin() const {
		return (const T*)data.data();
	}
	inline T* end() {
		return (T*)data.data() + right;
	}
	inline const T* end() const {
		return (const T*)data.data() + right;
	}
	inline T& front() {
		ASSERT(right > 0, "no data.");
		return data[0];
	}
	const inline T& front() const {
		ASSERT(right > 0, "no data.");
		return data[0];
	}
	inline T& back() {
		ASSERT(right > 0, "no data.");
		return data[right - 1];
	}
	const inline T& back() const {
		ASSERT(right > 0, "no data.");
		return data[right - 1];
	}
	inline bool contains(const T& value) const {
		for (const auto& dat : *this) {
			if (value == dat) return true;
		}
		return false;
	}
	inline vector<T> ToVector() {
		return vector<T>(begin(), end());
	}
	inline void Print() const {
		cout << '{';
		for (int i = 0; i < right; i++) cout << data[i] << ",}"[i == right - 1];
		cout << endl;
	}
	template<class S> inline auto AsType() const {
		return Stack<S, max_size>(*this);
	}
};

// ハッシュテーブル  // うまく実装できん
template<class T, int size = 0x100000, class KeyType = unsigned long long>
struct HashMap {
	array<pair<KeyType, T>, size> data;
	constexpr static KeyType mask = size - 1;
	constexpr static KeyType EMPTY = 0;
	constexpr static KeyType DELETED = 1;
	inline HashMap() {
		static_assert((size & size - 1) == 0, "not pow of 2");
		memset(&data[0], 0, sizeof(data));
	}
	const T& Get(const KeyType& key) const {
		// 既に値が格納されていることを仮定
		if (key == EMPTY || key == DELETED) return Get(key + (KeyType)2);
		auto address = key & mask;
		while (data[address].first != key) address = (address + 1) & mask;
		return data[address].second;
	}
	const T& Get(const KeyType& key, const T& default_value) const {
		// まだ値が格納されていない場合の値を指定
		if (key == EMPTY || key == DELETED) return Get(key + (KeyType)2, default_value);
		auto address = key & mask;
		while (true) {
			if (data[address].first == key) return data[address].second;
			else if (data[address].first == EMPTY) return default_value;
			address = (address + 1) & mask;
		}
	}
	void Set(const KeyType& key, const T& value) {
		// まだ値が格納されていないことを仮定
		if (key == EMPTY || key == DELETED) return Set(key + (KeyType)2, value);
		auto address = key & mask;
		while (data[address].first != EMPTY && data[address].first != DELETED) address = (address + 1) & mask;
		data[address].first = key;
		data[address].second = value;
	}
	T& operator[](const KeyType& key) {
		// 存在すればその値を返す
		// 存在しなければ新しく作って返す
		if (key == EMPTY || key == DELETED) return operator[](key + (KeyType)2);
		auto address = key & mask;
		while (data[address].first != EMPTY) {
			if (data[address].first == key) return data[address].second;
			address = (address + 1) & mask;
		}
		address = key & mask;
		while (data[address].first != EMPTY && data[address].first != DELETED) address = (address + 1) & mask;
		data[address].first = key;
		new(&data[address].second) T();
		return data[address].second;
	}
	void erase(const KeyType& key) {
		// 既に値が格納されていることを仮定
		if (key == EMPTY || key == DELETED) return erase(key + (KeyType)2);
		auto address = key & mask;
		while (data[address].first != key) address = (address + 1) & mask;
		data[address].first = DELETED;
		data[address].second.~T();  // これどうすればいいんだ
		memset(&data[address].second, 0, sizeof(T));
	}
	void clear() {
		// コストがでかい、デストラクタとか呼ばれない
		memset(&data[0], 0, sizeof(data));
	}

};
auto a = HashMap<double>();

// 時間 (秒)
inline double Time() {
	return static_cast<double>(chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now().time_since_epoch()).count()) * 1e-9;
}


// 重複除去
template<class VectorLike> inline void Deduplicate(VectorLike& vec) {
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
}


// 2 分法
template<class VectorLike, typename T> inline int SearchSorted(const VectorLike& vec, const T& a) {
	return lower_bound(vec.begin(), vec.end(), a) - vec.begin();
}


// argsort
template<typename T, int n, typename result_type, bool reverse = false> inline auto Argsort(const array<T, n>& vec) {
	array<result_type, n> res;
	iota(res.begin(), res.end(), 0);
	sort(res.begin(), res.end(), [&](const result_type& l, const result_type& r) {
		return reverse ? vec[l] > vec[r] : vec[l] < vec[r];
	});
	return res;
}


// popcount  // SSE 4.2 を使うべき
inline int Popcount(const unsigned int& x) {
#ifdef _MSC_VER
	return (int)__popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}
inline int Popcount(const unsigned long long& x) {
#ifdef _MSC_VER
	return (int)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
}

// x >> n & 1 が 1 になる最小の n ( x==0 は未定義 )
inline int CountRightZero(const unsigned int& x) {
#ifdef _MSC_VER
	unsigned long r;
	_BitScanForward(&r, x);
	return (int)r;
#else
	return __builtin_ctz(x);
#endif
}
inline int CountRightZero(const unsigned long long& x) {
#ifdef _MSC_VER
	unsigned long r;
	_BitScanForward64(&r, x);
	return (int)r;
#else
	return __builtin_ctzll(x);
#endif
}

inline double MonotonicallyIncreasingFunction(const double& h, const double& x) {
	// 0 < h < 1
	// f(0) = 0, f(1) = 1, f(0.5) = h
	ASSERT(h > 0.0 && h < 1.0, "0 < h < 1 not satisfied");
	if (h == 0.5) return x;
	const double& a = (1.0 - 2.0 * h) / (h * h);
	return expm1(log1p(a) * x) / a;
}
inline double MonotonicFunction(const double& start, const double& end, const double& h, const double& x) {
	// h: x = 0.5 での進捗率
	return MonotonicallyIncreasingFunction(h, x) * (end - start) + start;
}

#endif  // NAGISS_LIBRARY_HPP

// パラメータ
// K: 大きいほど未来の価値が小さくなる log2/100 = 0.007 くらいのとき野菜のインフレと釣り合う？

constexpr double K_START = 0.023684651315285833;  // OPTIMIZE [0.02, 0.06] LOG
constexpr double K_END = 0.04652107876916136;   // OPTIMIZE [0.01, 0.05] LOG
constexpr double K_H = 0.6676996005319429;      // OPTIMIZE [0.001, 0.999]

constexpr int hash_table_size = 19;
constexpr int beam_width = 200;

constexpr short PURCHASE_TURN_LIMIT = 828;  // OPTIMIZE [780, 880]

// 0 で通常
constexpr int SUBSCORE3_TIGHT_TURN = 0;     // OPTIMIZE [0, 2]

using ull = unsigned long long;
using i8 = int8_t;
using u8 = uint8_t;
using u16 = uint16_t;

// 入力
constexpr int N = 16;
constexpr int M = 5000;
constexpr int T = 1000;
array<i8, M> R;
array<i8, M> C;
array<u8, M> RC;
array<short, M> S;
array<short, M> E;
array<short, M> V;


struct alignas(32) BitBoard {
	__m256i data;
	union U {
		__m256i raveled;
		array<u16, 16> rows;
		array<ull, 4> ulls;
	};

	inline auto& Rows() {
		return ((U*)&data)->rows;
	}
	inline const auto& Rows() const {
		return ((U*)&data)->rows;
	}
	inline auto NonzeroIndices() const {
		auto res = Stack<u8, 256>();
		alignas(32) auto d = data;
		rep(i, 4) {
			auto& di = ((U*)&d)->ulls[i];
			while (di != 0ull) {
				const auto ctz = CountRightZero(di);
				res.push((u8)i * 64u + (u8)ctz);
				di ^= 1ull << ctz;
			}
		}
		return res;
	}
	inline bool Get(const u8& idx) const {
		return ((U*)&data)->ulls[idx >> 6] >> (idx & 63u) & 1u;
	}
	inline void Flip(const u8& idx) {
		((U*)&data)->ulls[idx >> 6] ^= 1ull << (idx & 63u);
	}
	inline BitBoard Left() const {
		return BitBoard{ _mm256_srli_epi16(data, 1) };
	}
	inline BitBoard Right() const {
		return BitBoard{ _mm256_slli_epi16(data, 1) };
	}
	inline BitBoard Up() const {
		return BitBoard{ _mm256_alignr_epi8(
			_mm256_permute2x128_si256(data, data, 0b10000001), data, 2  // _mm256_permute2x128_si256(data, data, 0b10000001) で data の上位ビットを取得
		) };  // alignr(s1, s2, mask) := ((s1 の上位と s2 の上位) >> mask) << 16 | ((s1 の下位と s2 の下位) >> mask) シフトは 8 bit 単位
	}
	inline BitBoard Down() const {
		return BitBoard{ _mm256_alignr_epi8(
			data, _mm256_permute2x128_si256(data, data, 0b00001000), 16 - 2  // _mm256_permute2x128_si256(data, data, 0b00001000) で data の下位ビットを上位に持ってくる
		) };
	}
	inline BitBoard& operator&=(const BitBoard& rhs) {
		data = _mm256_and_si256(data, rhs.data);
		return *this;
	}
	inline BitBoard& operator|=(const BitBoard& rhs) {
		data = _mm256_or_si256(data, rhs.data);
		return *this;
	}
	inline BitBoard& operator^=(const BitBoard& rhs) {
		data = _mm256_xor_si256(data, rhs.data);
		return *this;
	}
	inline BitBoard operator&(const BitBoard& rhs) const {
		return BitBoard(*this) &= rhs;
	}
	inline BitBoard operator|(const BitBoard& rhs) const {
		return BitBoard(*this) |= rhs;
	}
	inline BitBoard operator^(const BitBoard& rhs) const {
		return BitBoard(*this) ^= rhs;
	}
	inline BitBoard operator~() const {
		return BitBoard{ _mm256_xor_si256(data, _mm256_set1_epi32(-1)) };
	}
	inline bool Empty() const {
		return _mm256_testz_si256(data, data);
	}
	inline BitBoard& Expand() {  // 上下左右に広がる
		return *this |= Down() |= Right() |= Up() |= Left();
	}
	void Print() const {
		for (const auto& row : Rows()) {
			rep(i, 16) cout << (row >> i & 1) << ",\n"[i == 15];
		}
	}
};

namespace test {
void TestBitBoard() {
	auto bb = BitBoard{ 0 };
	bb.Flip(0 * 16 + 2);
	bb.Flip(2 * 16 + 4);
	bb.Flip(6 * 16 + 8);
	bb.Flip(255);
	bb.Print();

	cout << "empty" << endl;
	cout << bb.Empty() << endl;

	cout << "nonzero" << endl;
	bb.NonzeroIndices().AsType<short>().Print();

	cout << "expand" << endl;
	bb.Expand();
	bb.Print();

	cout << "left" << endl;
	bb = bb.Left();
	bb.Print();

	cout << "up" << endl;
	bb = bb.Up();
	bb.Print();

	cout << "reset" << endl;
	for (const auto& idx : bb.NonzeroIndices()) bb.Flip(idx);
	cout << "empty" << endl;
	cout << bb.Empty() << endl;
}
}


namespace board_index_functions {
u8 UpOf(u8 idx) {
	if (idx < 16u) return idx;
	else return idx - 16u;
}
u8 DownOf(u8 idx) {
	if (idx >= 240u) return idx;
	else return idx + 16u;
}
u8 RightOf(u8 idx) {
	if ((idx & 15u) == 15u) return idx;
	else return idx + 1u;
}
u8 LeftOf(u8 idx) {
	if ((idx & 15u) == 0u) return idx;
	else return idx - 1u;
}
}

namespace globals {
auto rng = Random(42);                              // random number generator
auto RNT = array<ull, 10000>();                     // random number table
auto EXP_NEG_KT = array<double, 1000>();
auto v_modified = array<double, M>();               // ターンで補正した野菜の価値
auto NEIGHBOR = array<array<BitBoard, 16>, 256>();  // neighbor[idx][d] := idx から距離 d 以内の場所たち
auto s_begins = array<short, T + 1>();              // t 日目の野菜の最初のインデックス
auto order_e = array<short, M>();                   // argsort(E)
auto e_begins = array<short, T + 1>();              // order_e のインデックスで、t 日目に消滅する野菜の最初のインデックス
auto start_bitboards = array<BitBoard, T>();        // そのターンに出現する野菜の位置
auto end_bitboards = array<BitBoard, T>();          // そのターンに消滅する野菜の位置
auto next_vegetable = array<short, M>();            // 同じマスに次に現れる野菜のインデックス  // 次がなければ -1

// ビームサーチ中に変動
auto t = 0;
auto future_value_table = Board<double, N, N>();   // 将来生える野菜の価値
auto current_money_table = Board<short, N, N>();   // 今生えてる野菜の価値 (補正なし)
auto current_index_table = Board<short, N, N>();   // 野菜インデックス
auto current_value_table = Board<double, N, N>();  // 今と将来の野菜の価値 (補正無し)
auto high_value_indices = array<u8, 256>();        // 今と将来の野菜の価値 (補正無し) をソートしたもの
auto next_end_table = Board<short, N, N>();        // 次にそのマスの価値が落ちるタイミング  // 次がなければ -1

void UpdateValueTable() {
	// State::Do をする前に呼ぶ
	
	// 出現
	rep3(idx_RCSEV, globals::s_begins[t], globals::s_begins[t + 1]) {
		const auto& rc = RC[idx_RCSEV];
		const auto& vm = v_modified[idx_RCSEV];
		const auto& v = V[idx_RCSEV];
		
		ASSERT(t == S[idx_RCSEV], "turn がおかしいよ");
		ASSERT(current_index_table.data[rc] < 0, "既に野菜があるよ");
		ASSERT(current_money_table.data[rc] == 0, "既に野菜があるよ");

		current_index_table.data[rc] = idx_RCSEV;
		current_money_table.data[rc] = v;
		future_value_table.data[rc] -= vm;

		ASSERT(future_value_table.data[rc] >= -1e3, "将来の価値がマイナスになることはないはずだよ");
	}

	// 消滅
	rep3(idx_order, e_begins[t], e_begins[t + 1]) {
		const auto& idx_RCSEV = order_e[idx_order];

		const auto& rc = RC[idx_RCSEV];
		const auto& v = V[idx_RCSEV];

		ASSERT(t == E[idx_RCSEV], "turn がおかしいよ");
		ASSERT(current_index_table.data[rc] == idx_RCSEV, "消滅させる野菜がないよ");
		ASSERT(current_money_table.data[rc] == v, "消滅させる野菜がないよ");

		current_index_table.data[rc] = -1;
		current_money_table.data[rc] = 0;

		ASSERT(next_end_table.data[rc] == t, "終わる turn 間違ってませんか");
		const auto& idx_next_vege = next_vegetable[idx_RCSEV];
		if (idx_next_vege != -1) {
			ASSERT(rc == RC[idx_next_vege], "場所間違ってませんか");
			next_end_table.data[rc] = E[idx_next_vege];
		}
		else {
			next_end_table.data[rc] = -1;
		}
	}

	t++;


	// これ 1 個先読みしないといけない…？？？？？うわ～～～～～～～
	// 1 turn 前ので近似… ごまかす…
	if (t != T) {
		rep(idx, 256) {
			current_value_table.data[idx] = current_money_table.data[idx] + future_value_table.data[idx] / EXP_NEG_KT[t];  // machine の数…は大丈夫だった
		}
		high_value_indices = Argsort<double, 256, u8, true>(current_value_table.data);
	}

}

}


struct State {
	BitBoard vegetables;
	BitBoard machines;  // 一定ターン以降は木を保つ
	short turn;  // 何も置いてない状態が 0
	short n_machines;
	int money;
	double score;
	double subscore2;
	double subscore3;

	unsigned hash;  // machines のみによって一意に定まる

	struct Action {
		// 新しく置くときは before == after にする
		u8 before, after;
	};
	struct NewStateInfo {
		double score;
		ull hash;
		Action action;
	};
	void Print() {
		cout << "State{" << endl;
		cout << "vegetables:" << endl;
		vegetables.Print();
		cout << "machines:" << endl;
		machines.Print();
		cout << "turn=" << turn << endl;
		cout << "n_machines=" << n_machines << endl;
		cout << "money=" << n_machines << endl;
		cout << "score=" << score << endl;
		cout << "subscore2=" << subscore2 << endl;
		cout << "subscore3=" << subscore3 << endl;
		cout << "}" << endl;
	}
	inline bool Terminated() const {
		return turn == 1000;
	}
	inline void Do(const Action& action) {
		// 1. 収穫機を移動させる
		//   - machine を変更する
		//   - future_value_table に応じて subscore2 を差分計算する
		//   - 野菜があれば money を増やす
		// 2. その日出現する野菜に応じて vegetables のビットを立てる
		// 3. machines と vegetables の共通部分を取って、野菜を収穫する
		//   - 収穫した部分の vegetables のビットは折る
		// 4. 野菜の出現に応じて future_value_table を減少させ、そのマスに machine があれば subscore2 を減らして money を増やす
		//   - 実際には変動があった場所のリストを作り、table は触らない
		// 5. その日消滅する野菜に応じて vegetables のビットを折る
		// 6. subscore3 を計算する
		//   - machine の無い、最も価値が高いマスについて、
		//     その価値が下がる前にそこに到達可能であれば、点をつける
		
		// 参照する外部の変数:
		// future_value_table
		// current_money_table
		// start_bitboards
		// end_bitboards
		// ほか


		// Step 1
		const auto& before = action.before;
		const auto& after = action.after;
		if (before == after) {
			if (machines.Get(after)) {
				// パスする場合
				// 何もしない
			}
			else {
				// 新しく置く場合
				n_machines++;
				money -= (int)n_machines * (int)n_machines * (int)n_machines;
				ASSERT(money >= 0, "money < 0");

				machines.Flip(after);
				hash ^= globals::RNT[after];
				subscore2 += globals::future_value_table.data[after];
				money += vegetables.Get(after) * n_machines * globals::current_money_table.data[after];
			}
		}
		else {
			// 移動させる場合
			ASSERT(machines.Get(before), "移動元に機械がないよ");
			ASSERT(!machines.Get(after), "移動先に機械があるよ");
			machines.Flip(before);
			machines.Flip(after);
			hash ^= globals::RNT[before] ^ globals::RNT[after];
			subscore2 += globals::future_value_table.data[after] - globals::future_value_table.data[before];
			money += vegetables.Get(after) * n_machines * globals::current_money_table.data[after];
		}

		// Step 2: 出現
		vegetables |= globals::start_bitboards[turn];

		// Step 3: 収穫(1)
		//auto intersection = machines & vegetables;
		//vegetables ^= intersection;
		vegetables.data = _mm256_andnot_si256(machines.data, vegetables.data);

		// Step 4: 収穫(2)
		rep3(idx_vegetables, globals::s_begins[turn], globals::s_begins[turn + 1]) {
			const auto& idx = RC[idx_vegetables];
			const auto& vm = globals::v_modified[idx_vegetables];
			const auto& v = V[idx_vegetables];
			
			if (machines.Get(idx)) {
				subscore2 -= vm;
				ASSERT(subscore2 >= -1e3, "subscore2 < 0");

				money += n_machines * v;
			}
		}

		// Step 5: 消滅
		vegetables.data = _mm256_andnot_si256(globals::end_bitboards[turn].data, vegetables.data);


		// !!!!!ターン増やすよ!!!!!
		turn++;
		// !!!!!ターン増やすよ!!!!!


		// Step 6: subscore3 の計算
		subscore3 = 0.0;
		
		remove_reference<decltype(globals::high_value_indices[0])>::type high_value_idx;  // 価値の高い場所
		for (int i = 0;; i++) {
			high_value_idx = globals::high_value_indices[i];
			if (!machines.Get(high_value_idx)) break;
		}
		if (globals::next_end_table.data[high_value_idx] != -1) {
			const auto value_decline_turn = min(
				globals::next_end_table.data[high_value_idx] - turn + 1 - SUBSCORE3_TIGHT_TURN,  // 価値が落ちるまでに行動できる回数  // 同じであっても 1 回猶予がある
				15
			);
			//ASSERT_RANGE(value_decline_turn, 1, 16);  // 先読みしてないのでこれにひっかかる…
			if (!(globals::NEIGHBOR[high_value_idx][value_decline_turn] & machines).Empty()) {  // 到達可能であれば
				subscore3 = globals::current_value_table.data[high_value_idx];
			}
		}
		

		if (turn == T) {
			score = money;
		}
		else {
			score = money
				+ (subscore2 / globals::EXP_NEG_KT[turn] + subscore3) * n_machines
				+ (int)((n_machines * (n_machines + 1)) / 2) * (int)((n_machines * (n_machines + 1)) / 2);  // 3 乗和
		}
	}
	template<class Vector>
	inline void GetNextStates(Vector& res) const {

		// TODO: n_machines が少ないときの処理

		if (((int)n_machines + 1) * ((int)n_machines + 1) * ((int)n_machines + 1) > money || turn >= PURCHASE_TURN_LIMIT) {
			// 資金が足りない場合 (1 個取り除く) or 一定ターン以降
			if (n_machines == 1) {
				// 機械が 1 個のとき
				const auto p_remove = machines.NonzeroIndices()[0];
				rep(p_add, 256) {
					if (p_remove == p_add) continue;
					auto new_state = *this;
					new_state.Do(Action{ p_remove, (u8)p_add });
					res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, (u8)p_add} });
				}
			}
			else {
				// 機械が 2 個以上のとき
				for (const auto& p_remove : machines.NonzeroIndices()) {
					// 葉じゃなかったら飛ばす
					using namespace board_index_functions;
					int neighbor_cnt = 0;  // 隣接する数
					for (const auto& drul : { DownOf(p_remove), RightOf(p_remove), UpOf(p_remove), LeftOf(p_remove) }) {
						if (drul != p_remove) {
							neighbor_cnt += machines.Get(drul);
						}
					}
					if (neighbor_cnt >= 2) continue;  // p は葉ではない

					// 実際に減らす
					auto machines_removed = machines;
					machines_removed.Flip(p_remove);

					// 1 個足す方法を探す
					// 元々の連結成分に 1 箇所で隣接 <=> xor が 1 かつ, 2 ペアの or の xor が 1
					auto&& cand = (machines_removed.Down() ^ machines_removed.Right() ^ machines_removed.Up() ^ machines_removed.Left())
						& ((machines_removed.Down() | machines_removed.Right()) ^ (machines_removed.Up() | machines_removed.Left()))
						& ~machines;
					for (const auto& p_add : cand.NonzeroIndices()) {
						//if (p_remove == p_add) continue;
						ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");
						auto new_state = *this;
						new_state.Do(Action{ p_remove, p_add });
						res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
					}
				}
			}
		}
		else {
			// 資金が足りてる場合
			if (n_machines == 0) {
				rep(p_add, 256) {
					auto new_state = *this;
					new_state.Do(Action{ (u8)p_add, (u8)p_add });
					res.push(NewStateInfo{ new_state.score, new_state.hash, {(u8)p_add, (u8)p_add} });
				}
			}
			else {
				// 1 個足す方法を探す
				auto&& cand = (machines.Down() ^ machines.Right() ^ machines.Up() ^ machines.Left())
					& ((machines.Down() | machines.Right()) ^ (machines.Up() | machines.Left()))
					& ~machines;
				for (const auto& p_add : cand.NonzeroIndices()) {
					auto new_state = *this;
					new_state.Do(Action{ p_add, p_add });
					res.push(NewStateInfo{ new_state.score, new_state.hash, {p_add, p_add} });
				}
			}
		}

		// 角で v が w になるやつとかも考慮すべき？

	}

};

void Solve() {
	// 入力を受け取る
	{
		int buf;
		scanf("%d %d %d", &buf, &buf, &buf);
		rep(i, M) {
			scanf("%hhd %hhd %hd %hd %hd", &R[i], &C[i], &S[i], &E[i], &V[i]);
			RC[i] = ((u8)R[i] * (u8)N + (u8)C[i]);
		}
	}

	// 色々初期化
	{
		using namespace globals;
		for (auto&& r : RNT) {
			r = (unsigned)rng.next();
		}
		EXP_NEG_KT[0] = 1.0;
		rep1(i, T - 1) EXP_NEG_KT[i] = EXP_NEG_KT[i - 1] * exp(-MonotonicFunction(K_START, K_END, K_H, (double)i / (double)T));
		rep(i, M) {
			v_modified[i] = V[i] * EXP_NEG_KT[S[i]];
			future_value_table[{ R[i], C[i] }] += v_modified[i];
			start_bitboards[S[i]].Flip(RC[i]);
			end_bitboards[E[i]].Flip(RC[i]);
		}
		// next_vegetable
		auto next_vegetable_board = Board<short, N, N>();
		next_vegetable_board.Fill(-1);
		next_end_table.Fill(-1);
		for (int i = M - 1; i >= 0; i--) {
			next_vegetable[i] = next_vegetable_board.data[RC[i]];
			next_vegetable_board.data[RC[i]] = i;
			next_end_table.data[RC[i]] = E[i];
		}
		// NEIGHBOR
		rep(i, 256) {
			NEIGHBOR[i][0].Flip((u8)i);
			rep(d, (ll)NEIGHBOR[i].size() - 1) {
				NEIGHBOR[i][d + 1] = NEIGHBOR[i][d];
				NEIGHBOR[i][d + 1].Expand();
			}
		}

		// s_begins
		auto idx = (short)0;
		s_begins[0] = 0;
		rep(turn, T) {
			while (idx < (int)S.size() && S[idx] <= turn) {
				idx++;
			}
			s_begins[turn + 1] = idx;
		}
		// order_e
		order_e = Argsort<short, M, short>(E);
		// e_begins
		idx = (short)0;
		e_begins[0] = 0;
		rep(turn, T) {
			while (idx < (int)E.size() && E[order_e[idx]] <= turn) {
				idx++;
			}
			e_begins[turn + 1] = idx;
		}

		current_index_table.Fill(-1);

		high_value_indices = Argsort<double, 256, u8, true>(future_value_table.data);
	}

	//globals::future_value_table.Print();
	//cout << "NEIGHBOR[200][7]" << endl;
	//globals::NEIGHBOR[200][7].Print();

	// ビームサーチ
	
	{
		struct Node {
			double score;
			Node* parent_node;
			State* state;
			typename State::Action action;

			inline bool operator<(const Node& rhs) const { return score < rhs.score; }
			inline bool operator>(const Node& rhs) const { return score > rhs.score; }
		};
		static Stack<State, 1 + beam_width * T> state_buffer;
		static Stack<Node, 1 + beam_width * T> node_buffer;
		state_buffer.push(State{});
		state_buffer.back().money = 1;
		node_buffer.push({ state_buffer[0].score, nullptr, &state_buffer[0] });
		static Stack<Node, 200000> q;
		Node* parent_nodes_begin = node_buffer.begin();
		Node* parent_nodes_end = node_buffer.end();
		Node* best_node = nullptr;

		rep(t, T) {
			static HashMap<Node*, 1 << hash_table_size> dict_hash_to_candidate;
			dict_hash_to_candidate.clear();

			for (auto parent_node = parent_nodes_begin; parent_node != parent_nodes_end; parent_node++) {
				static Stack<typename State::NewStateInfo, 10000> next_states;
				next_states.clear();
				parent_node->state->GetNextStates(next_states);
				for (const auto& r : next_states) {
					auto& old_node = dict_hash_to_candidate[r.hash];
					if (old_node == NULL) {
						// まだそのハッシュの状態が無い
						q.push({ r.score, parent_node, nullptr, r.action });
						old_node = &q.back();
					} else if (old_node->score < r.score) {
						// 同じハッシュでスコアがより良い
						old_node->score = r.score;
						old_node->parent_node = parent_node;
						old_node->action = r.action;
					}
				}

			}
			cerr << "q.size()=" << q.size() << endl;

			if (beam_width < q.size()) {
				nth_element(q.begin(), q.begin() + beam_width, q.end(), greater<>());  // ここでハッシュテーブル破壊されている
				q.resize(beam_width);
			}

			for (auto&& node : q) {
				auto state = *node.parent_node->state;
				state.Do(node.action);
				state_buffer.push(state);
				node.state = &state_buffer.back();
				node_buffer.push(node);
				if (state.Terminated()) {
					if (best_node == nullptr || best_node->score < state.score) {
						best_node = &node_buffer.back();
					}
				}
			}
			q.clear();
			parent_nodes_begin = parent_nodes_end;
			parent_nodes_end = node_buffer.end();

			globals::UpdateValueTable();
		}

		// 結果を出力
		{
			ASSERT(best_node != nullptr, "best_node がないよ");
			auto path = array<State::Action, T>();
			auto node = best_node;
			rep(i, T) {
				path[T - 1 - i] = node->action;
				node = node->parent_node;
			}
			ASSERT(node->parent_node == nullptr, "根ノードじゃないよ");
			for (const auto& action : path) {
				const auto& before = action.before;
				const auto& after = action.after;
				if (before == after) {
					cout << (short)(after >> 4) << " " << (short)(after & 15u) << endl;
				}
				else {
					cout << (short)(before >> 4) << " " << (short)(before & 15u) << " "
						<< (short)(after >> 4) << " " << (short)(after & 15u) << endl;
				}
			}
			cerr << best_node->score << endl;
		}
	}

	
}


int main() {
	//test::TestBitBoard();
	Solve();

}

#ifdef __GNUC__
#pragma clang attribute pop
#endif


/*
- 終盤のインフレがすごいが終盤はあまり動けない
- 最重要: subscore3 の実装
- ビームサーチの時間調整
- ハッシュを雑に
- 価値の低い野菜の無視
- ハッシュが重いしいらないかもしれない
- ビーム幅 200 からの候補 60000, 多すぎる
*/

