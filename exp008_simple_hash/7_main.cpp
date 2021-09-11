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


template<class T, int size = 0x100000, class KeyType = unsigned>
struct MinimumHashMap {
	// ハッシュの値が size 以下
	array<T, size> data;
	Stack<int, size> used;
	constexpr static KeyType mask = size - 1;
	inline MinimumHashMap() {
		static_assert((size & size - 1) == 0, "not pow of 2");
		memset(&data[0], (unsigned char)-1, sizeof(data));
	}
	inline T& operator[](const KeyType& key) {
		if (data[key] == (T)-1) used.push(key);
		return data[key];
	}
	inline void clear() {
		for (const auto& key : used) data[key] = (T)-1;
		used.right = 0;
	}
};

// ハッシュテーブル  // うまく実装できん
template<class T, int size = 0x100000, class KeyType = unsigned long long>
struct SlowHashMap {
	array<pair<KeyType, T>, size> data;
	constexpr static KeyType mask = size - 1;
	constexpr static KeyType EMPTY = 0;
	constexpr static KeyType DELETED = 1;
	inline SlowHashMap() {
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


// P 制御
struct PController {
	double Kp;
	inline PController(const double& kp) : Kp(kp) {}
	inline double operator()(const double& deviation){
		return Kp * deviation;
	}
};

// PID 制御
struct PIDController {
	double Kp, Ki, Kd;
	double sum_deviation, old_deviation;
	inline PIDController(const double& kp, const double& ki, const double& kd)
		: Kp(kp), Ki(ki), Kd(kd), sum_deviation(0.0), old_deviation(0.0) {}
	inline double operator()(const double& deviation) {
		sum_deviation += deviation;
		const auto p = Kp * deviation;
		const auto i = Ki * sum_deviation;
		const auto d = Kd * (deviation - old_deviation);
		old_deviation = deviation;
		cerr << "pid: " << p << " " << i << " " << d << "\n";
		return p + i + d;
	}
};

#endif  // NAGISS_LIBRARY_HPP

// パラメータ

#ifdef _MSC_VER
constexpr double TIME_LIMIT = 4.0;
#else
constexpr double TIME_LIMIT = 1.7;
#endif
constexpr int hash_table_size = 9;         // OPTIMIZE [9, 18]


// K: 大きいほど未来の価値が小さくなる log2/100 = 0.007 くらいのとき野菜のインフレと釣り合う？
constexpr double K_START = 0.056766944944961706;  // OPTIMIZE [0.04, 0.2] LOG
constexpr double K_END = 0.05225273209787457;   // OPTIMIZE [0.03, 0.1] LOG
constexpr double K_H = 0.8085952713056656;      // OPTIMIZE [0.001, 0.999]

constexpr short PURCHASE_TURN_LIMIT = 833;  // OPTIMIZE [790, 870]

// 0 で通常
constexpr int SUBSCORE3_TIGHT_TURN = 0;     // OPTIMIZEd

constexpr int ROUGH_HASH = 0;      // OPTIMIZE {0, 0b00000001, 0b00010001, 0b00010011, 0b00110011}

// ビーム
constexpr double TARGET_BEAM_WIDTH_INCREASE_RATE = 3.566604338929269;      // OPTIMIZE [0.25, 4.0] LOG
constexpr double TARGET_BEAM_WIDTH_HALF_PROGRES_RATE = 0.5640714305297574;  // OPTIMIZE [0.02, 0.98]
constexpr auto MAX_BEAM_WIDTH = 682;                        // OPTIMIZE [400, 4000] LOG
constexpr auto MIN_BEAM_WIDTH = 50;

// 型
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
auto T0 = Time();
auto rng = Random(123456789);                       // random number generator
auto RNT = array<unsigned, 10000>();                     // random number table
auto EXP_NEG_KT = array<double, 1000>();
auto v_modified = array<double, M>();               // ターンで補正した野菜の価値
auto NEIGHBOR = array<array<BitBoard, 16>, 256>();  // neighbor[idx][d] := idx から距離 d 以内の場所たち
auto s_begins = array<short, T + 1>();              // t 日目の野菜の最初のインデックス
auto order_e = array<short, M>();                   // argsort(E)
auto e_begins = array<short, T + 1>();              // order_e のインデックスで、t 日目に消滅する野菜の最初のインデックス
auto start_bitboards = array<BitBoard, T>();        // そのターンに出現する野菜の位置
auto end_bitboards = array<BitBoard, T>();          // そのターンに消滅する野菜の位置
auto next_vegetable = array<short, M>();            // 同じマスに次に現れる野菜のインデックス  // 次がなければ -1
auto exact_future_value = array<double, M>();       // その野菜が出現した後のそのマスの future_value の値

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
		//future_value_table.data[rc] -= vm;
		future_value_table.data[rc] = exact_future_value[idx_RCSEV];

		ASSERT(future_value_table.data[rc] >= 0.0, "将来の価値がマイナスになることはないはずだよ");
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
		unsigned hash;
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
		// 7. subscore4 を計算する
		//   - 斜めが多いと将来のスコアが減る
		
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
				hash += globals::RNT[after | ROUGH_HASH];
				hash &= (1 << hash_table_size) - 1;
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
			hash += globals::RNT[after | ROUGH_HASH] - globals::RNT[before | ROUGH_HASH];
			hash &= (1 << hash_table_size) - 1;
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
		rep(i, 32) {
			high_value_idx = globals::high_value_indices[i];
			if (vegetables.Get(high_value_idx)) break;
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
		
		// Step 7
		auto penalty = 1.0;
		/*
		if (!(machines.Left() & BitBoard { _mm256_andnot_si256(machines.data, machines.Right().data) }).Empty()) {
			penalty = 0.98;
		}
		*/


		if (turn == T) {
			score = money;
		}
		else {
			score = money
				+ (subscore2 / globals::EXP_NEG_KT[turn] + subscore3) * n_machines * penalty
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
						//if (turn >= 50 && turn < T - 1
						//	&& new_state.score - new_state.subscore3 * new_state.n_machines < (score - subscore3 * n_machines) * 0.8) continue;  // 枝刈り  // 悪化
						res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
					}
				}

				// 動かす特殊ケース v -> w
				const auto l = machines.Left();
				const auto r = machines.Right();
				const auto u = machines.Up();
				const auto d = machines.Down();
				const auto dd = d.Down();
				const auto dl = d.Left();
				const auto ddl = dd.Left();
				const auto dll = dl.Left();
				const auto dr = d.Right();
				const auto ddr = dd.Right();
				const auto drr = dr.Right();
				const auto ur = u.Right();
				const auto rr = r.Right();
				{
					using namespace board_index_functions;
					const auto cond_center = d & ~u;
					{
						const auto cond_dl = cond_center & l & ~r & ~dll & ~ddl;
						// 左下から右上
						for (const auto& p_remove : (cond_dl& machines & ~dl).NonzeroIndices()) {
							const auto& p_add = UpOf(RightOf(p_remove));
							ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
						}
						// 右上から左下
						for (const auto& p_add : (cond_dl & ~machines & dl).NonzeroIndices()) {
							const auto& p_remove = UpOf(RightOf(p_add));
							ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
						}
					}
					{
						const auto cond_dr = cond_center & r & ~l & ~drr & ~ddr;
						// 右下から左上
						for (const auto& p_remove : (cond_dr & machines & ~dr).NonzeroIndices()) {
							const auto p_add = UpOf(LeftOf(p_remove));
							ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
						}
						// 右上から左下
						for (const auto& p_add : (cond_dr & ~machines & dr).NonzeroIndices()) {
							const auto p_remove = UpOf(LeftOf(p_add));
							ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
						}
					}
				}
				// 特殊ケース 2: [ -> ]
				{
					using namespace board_index_functions;
					const auto cond_ud = l & r & dl & dr & ~u & ~dd;
					// 上から下
					for (const auto& p_add : (cond_ud & ~machines & d).NonzeroIndices()) {
						const auto p_remove = UpOf(p_add);
						ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
					}
					// 下から上
					for (const auto& p_remove : (cond_ud & machines & ~d).NonzeroIndices()) {
						const auto p_add = UpOf(p_remove);
						ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
					}
					const auto cond_lr = u & d & ur & dr & ~l & ~rr;
					// 左から右
					for (const auto& p_add : (cond_lr & ~machines & r).NonzeroIndices()) {
						const auto p_remove = LeftOf(p_add);
						ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
					}
					// 右から左
					for (const auto& p_remove : (cond_lr & machines & ~r).NonzeroIndices()) {
						const auto p_add = LeftOf(p_remove);
						ASSERT(p_remove != p_add, "元と同じ箇所は選ばれないはずだよ");  auto new_state = *this;  new_state.Do(Action{ p_remove, p_add });  res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
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

namespace beam_width_control {

// 変動なし
constexpr auto BASE_SEC_PER_WIDTH = array<int, T>{14048,14035,14011,13974,13927,13869,13801,13725,13641,13551,13456,13358,13257,13156,13055,12955,12858,12765,12676,12593,12516,12445,12380,12323,12272,12229,12192,12161,12137,12119,12106,12098,12095,12095,12100,12108,12119,12132,12147,12164,12181,12200,12218,12236,12254,12271,12286,12299,12311,12320,12327,12331,12332,12330,12325,12318,12309,12297,12284,12269,12253,12238,12222,12208,12195,12184,12175,12170,12167,12168,12172,12180,12191,12204,12221,12240,12261,12283,12306,12330,12354,12377,12400,12422,12443,12463,12481,12499,12516,12532,12547,12563,12578,12593,12609,12625,12642,12659,12678,12697,12716,12736,12757,12778,12799,12821,12842,12863,12885,12906,12928,12949,12971,12993,13016,13039,13063,13088,13115,13142,13171,13201,13233,13266,13300,13335,13371,13407,13444,13480,13517,13553,13588,13622,13655,13687,13717,13746,13773,13798,13823,13846,13867,13888,13908,13928,13947,13965,13984,14002,14021,14039,14058,14076,14095,14114,14132,14150,14168,14185,14202,14218,14234,14249,14263,14276,14289,14301,14312,14323,14334,14345,14356,14367,14378,14390,14403,14417,14432,14448,14465,14484,14504,14525,14547,14571,14595,14620,14646,14673,14700,14727,14755,14782,14810,14837,14864,14890,14916,14941,14966,14990,15014,15036,15059,15080,15102,15123,15144,15165,15186,15208,15230,15252,15276,15300,15325,15350,15377,15404,15432,15461,15490,15519,15549,15578,15606,15635,15662,15688,15714,15738,15761,15783,15804,15824,15843,15861,15878,15895,15911,15927,15942,15957,15971,15985,15999,16012,16025,16037,16048,16060,16071,16081,16092,16103,16115,16127,16140,16155,16172,16190,16211,16234,16259,16287,16318,16350,16386,16423,16462,16502,16544,16586,16628,16670,16712,16753,16792,16830,16866,16901,16933,16962,16990,17014,17037,17057,17075,17091,17105,17117,17128,17137,17144,17151,17156,17161,17166,17170,17174,17177,17181,17186,17190,17196,17201,17207,17213,17220,17227,17234,17241,17247,17253,17258,17263,17266,17269,17270,17270,17269,17267,17263,17259,17254,17250,17245,17240,17237,17235,17235,17237,17243,17251,17263,17279,17300,17325,17355,17390,17429,17473,17523,17576,17635,17697,17764,17834,17907,17984,18063,18144,18227,18311,18396,18482,18568,18653,18738,18822,18904,18985,19064,19141,19216,19288,19358,19425,19489,19550,19608,19663,19713,19760,19803,19842,19876,19906,19931,19952,19967,19978,19984,19986,19983,19977,19967,19955,19939,19922,19903,19883,19862,19842,19823,19805,19789,19775,19764,19756,19751,19751,19755,19763,19775,19793,19816,19844,19877,19915,19959,20008,20063,20122,20188,20258,20333,20413,20497,20585,20676,20770,20866,20964,21061,21158,21254,21348,21438,21524,21606,21683,21753,21818,21876,21928,21974,22014,22049,22078,22103,22124,22141,22157,22170,22182,22193,22204,22215,22227,22239,22252,22266,22281,22296,22312,22328,22344,22361,22377,22394,22410,22426,22443,22459,22476,22494,22513,22533,22556,22582,22611,22644,22681,22723,22770,22822,22880,22943,23012,23086,23165,23248,23335,23425,23518,23612,23708,23803,23899,23994,24086,24177,24265,24350,24432,24511,24585,24656,24724,24787,24848,24905,24960,25013,25063,25112,25161,25209,25258,25308,25359,25413,25469,25529,25593,25661,25733,25810,25891,25977,26068,26162,26259,26360,26463,26567,26672,26777,26880,26982,27082,27178,27271,27359,27443,27522,27595,27664,27728,27786,27841,27891,27938,27982,28022,28061,28098,28133,28167,28201,28234,28266,28299,28332,28364,28396,28429,28461,28492,28523,28554,28584,28614,28644,28673,28703,28734,28765,28799,28835,28874,28918,28966,29020,29080,29147,29222,29304,29395,29494,29601,29715,29838,29967,30102,30242,30387,30536,30687,30839,30992,31145,31296,31446,31592,31736,31875,32011,32141,32267,32388,32503,32612,32717,32816,32910,32998,33082,33162,33237,33309,33378,33445,33509,33572,33634,33696,33758,33821,33885,33951,34019,34089,34161,34237,34315,34395,34479,34566,34655,34748,34843,34941,35042,35145,35252,35361,35473,35587,35703,35822,35942,36063,36185,36307,36429,36550,36671,36790,36909,37025,37140,37254,37367,37479,37591,37703,37817,37931,38048,38168,38290,38415,38544,38676,38810,38948,39088,39229,39372,39515,39658,39800,39941,40080,40216,40351,40482,40611,40737,40861,40982,41101,41217,41332,41444,41554,41662,41768,41871,41972,42070,42166,42259,42350,42438,42523,42607,42690,42772,42853,42935,43017,43102,43189,43279,43374,43473,43577,43686,43801,43921,44048,44179,44315,44456,44601,44749,44899,45052,45205,45359,45513,45666,45819,45971,46123,46274,46426,46578,46731,46887,47046,47207,47373,47544,47720,47900,48086,48277,48473,48672,48875,49081,49289,49497,49705,49912,50118,50320,50519,50714,50905,51092,51274,51452,51625,51794,51960,52122,52281,52438,52593,52746,52898,53049,53200,53350,53500,53651,53802,53955,54108,54263,54419,54578,54738,54901,55065,55232,55402,55573,55746,55921,56097,56274,56451,56628,56803,56977,57148,57317,57481,57642,57798,57950,58097,58240,58380,58517,58651,58786,58920,59057,59198,59345,59499,59663,59839,60027,60231,60452,60692,60951,61232,61536,61862,62212,62585,62981,63399,63837,64292,64763,65246,65737,66232,66726,67215,67693,68156,68597,69014,69401,69755,70074,70355,70598,70802,70969,71100,71199,71268,71313,71337,71346,71344,71336,71327,71319,71317,71323,71339,71364,71400,71445,71498,71558,71622,71689,71755,71819,71879,71933,71981,72020,72051,72074,72090,72099,72102,72101,72097,72092,72085,72080,72076,72075,72076,72080,72088,72098,72110,72124,72139,72155,72170,72184,72196,72207,72215,72221,72225,72227,72227,72227,72225,72224,72223,72224,72226,72230,72237,72246,72258,72272,72288,72306,72325,72345,72364,72383,72401,72416,72428,72438,72443,72445,72443,72436,72426,72412,72394,72374,72352,72329,72305,72282,72261,72242,72227,72216,72208,72206,72208,72215,72225,72238,72253,72269,72285,72298,72309,72314,72314,72307,72291,72266,72231,72186,72131,72064,71986,71897,71796,71683,71559,71423,71275,71115,70943,70759,70562,70354,70134,69902,69659,69406,69143,68871,68591,68304,68010,67712,67409,67103,66795,66485,66175,65866,65558,65253,64952,64657,64368,64088,63818,63561,63317,63090,62881,62693,62527,62386,62270,62183,62123,62094};
double beam_search_time_limit;
double t_beam_search;
double mean_expected_base_sec;
int TURN_FIX = 950;


// 変動あり
int turn = 0;
double cum_base_sec = 0.0;  // s
double remaining_expected_base_sec = 0.0;  // b
int beam_width_at_turn_fix;

inline double Schedule(const double& t) {
	// パラメータ
	// 相対的なビーム幅の変化
	return MonotonicFunction(1.0, TARGET_BEAM_WIDTH_INCREASE_RATE, TARGET_BEAM_WIDTH_HALF_PROGRES_RATE, t);
	//return 1.0;
}
inline double ExpectedBaseSec(const int& turn_) {  // c
	return (double)BASE_SEC_PER_WIDTH[turn_] * Schedule((double)turn_ / (double)T);
}
int BeamWidth() {
	// ループの最初で呼ぶ
	// 呼ばれるたびに状態を更新する
	if (turn == 0) {
		// 初期化
		for (int i = 0; i < T; i++) remaining_expected_base_sec += ExpectedBaseSec(i);
		mean_expected_base_sec = remaining_expected_base_sec / (double)T;
		t_beam_search = Time();
		beam_search_time_limit = TIME_LIMIT - (t_beam_search - globals::T0);
	}

	// 幅を求める
	int beam_width;
	{
		const auto elapsed_time = Time() - t_beam_search;  // t
		const auto remaining_time = beam_search_time_limit - elapsed_time;  // r

		// 最初用の処理
		const auto additional_elapsed_time = beam_search_time_limit * 0.005;
		const auto additional_cum_base_sec = 1e9 * additional_elapsed_time;

		// 最後用の処理
		/*
		const auto additional_remaining_time = 0.2;  // r'
		//const auto additional_remaining_expected_base_sec = additional_remaining_time * (cum_base_sec + additional_cum_base_sec) / (elapsed_time + additional_elapsed_time);  // これだめだ！！！

		const auto base_sec_processing_time = (elapsed_time + additional_elapsed_time) / (cum_base_sec + additional_cum_base_sec);  // e
		const auto remaining_processable_base_sec = remaining_time / base_sec_processing_time;                         // v
		const auto additional_remaining_processable_base_sec = additional_remaining_time / base_sec_processing_time;   // v'

		// 
		const auto base_sec = (remaining_processable_base_sec + additional_remaining_processable_base_sec)
			                * ExpectedBaseSec(turn) / (remaining_expected_base_sec + mean_expected_base_sec * additional_remaining_time * T / beam_search_time_limit);  // q = v * c / b
		beam_width = base_sec / (double)BASE_SEC_PER_WIDTH[turn];

		//beam_width = clipped()

		// src / tbu
		//beam_width = ((cum_base_sec + additional_cum_base_sec) * (remaining_time              + additional_remaining_time)              * ExpectedBaseSec(turn))
		//	       / ((elapsed_time + additional_elapsed_time) * (remaining_expected_base_sec + additional_remaining_expected_base_sec) * (double)BASE_SEC_PER_WIDTH[turn]);

		cerr << "elapsed_time=" << elapsed_time << "  cum_base_sec/elapsed_time=" << cum_base_sec / elapsed_time << "  modified=" << (cum_base_sec + additional_cum_base_sec) / (elapsed_time + additional_elapsed_time) << "\n";
		cerr << "remaining_time=" << remaining_time << "\n";
		// これうそ！！
		//cerr << "remaining_expected_base_sec=" << remaining_expected_base_sec << "  expected_remaining_time=" << remaining_expected_base_sec / (cum_base_sec / elapsed_time) << "\n";
		*/

		// 綺麗な方法わからん！！！！！！
		if (turn <= TURN_FIX) {
			const auto base_sec_processing_time = (elapsed_time + additional_elapsed_time) / (cum_base_sec + additional_cum_base_sec);  // e  // 1e-9 くらい
			const auto remaining_processable_base_sec = remaining_time / base_sec_processing_time;                         // v = r / e
			const auto base_sec = remaining_processable_base_sec * ExpectedBaseSec(turn) / remaining_expected_base_sec;  // q = v * c / b
			beam_width = base_sec / (double)BASE_SEC_PER_WIDTH[turn];  // q / u
			if (turn == TURN_FIX) {
				beam_width_at_turn_fix = beam_width;
			}
		}
		else {
			beam_width = beam_width_at_turn_fix * Schedule(turn) / Schedule(TURN_FIX);
		}
		beam_width = clipped(beam_width, MIN_BEAM_WIDTH, MAX_BEAM_WIDTH);
		if (remaining_time < -0.1) beam_width = 24;
		if (turn % 50 == 49) {
			cerr << "elapsed_time=" << elapsed_time << "  cum_base_sec/elapsed_time=" << cum_base_sec / elapsed_time << "  modified=" << (cum_base_sec + additional_cum_base_sec) / (elapsed_time + additional_elapsed_time) << "\n";
			cerr << "remaining_time=" << remaining_time << "\n";
		}
	}
	cum_base_sec += (double)beam_width * (double)BASE_SEC_PER_WIDTH[turn];
	remaining_expected_base_sec -= ExpectedBaseSec(turn);
	turn++;
	return beam_width;
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
			r = (unsigned)rng.next() & (1 << hash_table_size) - 1;
		}
		EXP_NEG_KT[0] = 1.0;
		rep1(i, T - 1) EXP_NEG_KT[i] = EXP_NEG_KT[i - 1] * exp(-MonotonicFunction(K_START, K_END, K_H, (double)i / (double)T));
		rep(i, M) {
			v_modified[i] = V[i] * EXP_NEG_KT[S[i]];
			start_bitboards[S[i]].Flip(RC[i]);
			end_bitboards[E[i]].Flip(RC[i]);
		}
		// future_value_table, next_vegetable
		auto next_vegetable_board = Board<short, N, N>();
		next_vegetable_board.Fill(-1);
		next_end_table.Fill(-1);
		for (int i = M - 1; i >= 0; i--) {
			exact_future_value[i] = future_value_table.data[RC[i]];
			future_value_table.data[RC[i]] += v_modified[i];
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
		
		//auto beam_width_controller = PController(500.0);
		//auto beam_width_controller = PIDController(100.0, 1.0, 1000.0);
		struct Node {
			double score;
			Node* parent_node;
			State* state;
			typename State::Action action;

			inline bool operator<(const Node& rhs) const { return score < rhs.score; }
			inline bool operator>(const Node& rhs) const { return score > rhs.score; }
		};

		static Stack<State, 1 + MAX_BEAM_WIDTH * T> state_buffer;
		static Stack<Node, 1 + MAX_BEAM_WIDTH * T> node_buffer;
		state_buffer.push(State{});
		state_buffer.back().money = 1;
		node_buffer.push({ state_buffer[0].score, nullptr, &state_buffer[0] });
		static Stack<Node, 500000> q;
		Node* parent_nodes_begin = node_buffer.begin();
		Node* parent_nodes_end = node_buffer.end();
		Node* best_node = nullptr;
		const double t_beam_search = Time();
		//const double beam_search_time_limit = TIME_LIMIT - (t_beam_search - globals::T0);
		int beam_width;
		rep(t, T) {
			static MinimumHashMap<int, 1 << hash_table_size> dict_hash_to_candidate;
			dict_hash_to_candidate.clear();

			beam_width = beam_width_control::BeamWidth();
			for (auto parent_node = parent_nodes_begin; parent_node != parent_nodes_end; parent_node++) {
				static Stack<typename State::NewStateInfo, 10000> next_states;
				next_states.clear();
				parent_node->state->GetNextStates(next_states);
				for (const auto& r : next_states) {
					//if (t >= 50 && t < T - 1 && r.score < parent_node->score * 0.9) continue;  // 枝刈り  // 悪化
					auto& idx_old_node = dict_hash_to_candidate[r.hash];
					if (idx_old_node == -1) {
						// まだそのハッシュの状態が無い
						idx_old_node = q.size();
						q.push({ r.score, parent_node, nullptr, r.action });
					} else if (q[idx_old_node].score < r.score) {
						// 同じハッシュでスコアがより良い
						q[idx_old_node].score = r.score;
						q[idx_old_node].parent_node = parent_node;
						q[idx_old_node].action = r.action;
					}
				}

			}


			// ビーム幅制御
			{
				/* 没
				const auto real_progress = (double)t / (double)T;
				const auto target_progress = (Time() - t_beam_search) / beam_search_time_limit;
				beam_width *= (200 + beam_width_controller(real_progress - target_progress)) / 200;
				beam_width = clipped(beam_width, MIN_BEAM_WIDTH, MAX_BEAM_WIDTH);
				*/

			}
			if (true) {
				static double t_last = t_beam_search;

				if (t % 50 == 49) {
					//cerr << "real_progress=" << real_progress << "\n";
					//cerr << "target_progress=" << target_progress << "\n";
					//cerr << "deviation=" << real_progress - target_progress << "\n";  // 大きい -> 速く進めすぎ -> ビーム幅を大きくする
					cerr << "turn=" << t << "\n";
					//cerr << "time=" << Time() - t_beam_search << "\n";
					cerr << "time / width = " << (Time() - t_last) / beam_width << endl;
					cerr << "beam_width=" << beam_width << "\n";
					cerr << "q.size()=" << q.size() << "\n";
					cerr << endl;
				}
				t_last = Time();
			}

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

			// ~~~~~~~ 統計用に幅あたりの処理時間を出力する ~~~~~~~
			{
				//static double t_last = t_beam_search;
				//cerr << (Time() - t_last) / beam_width << "\n";
				//t_last = Time();
			}
			// ~~~~~~~
		}

		// 結果を出力
		cerr << best_node->score << endl;
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
- 最重要: ビームサーチの時間調整
- 重要: ハッシュを雑に
  - 微妙か？評価関数改善して誘導したほうが強そう
- 価値の低い野菜の無視
- ハッシュが重いしいらないかもしれない
- ビーム幅 200 からの候補 60000, 多すぎる
- subscore3 の改善
- 前 turn より減ってたら採用しない感じの枝刈り
- 斜めが少ないほど良い？
*/

