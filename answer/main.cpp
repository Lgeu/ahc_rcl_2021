#ifndef NAGISS_LIBRARY_HPP
#define NAGISS_LIBRARY_HPP
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
#include<cstdlib>
#include<ctime>
#include<string>
#include<sstream>
#include<chrono>
#include<climits>
#include<intrin.h>

#ifdef __GNUC__
//#pragma GCC target("avx2")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native")
#pragma GCC optimize("O3")
#pragma GCC optimize("Ofast")
//#pragma GCC optimize("unroll-loops")
#include<x86intrin.h>
#endif

// ========================== macroes ==========================

#define rep(i,n) for(ll (i)=0; (i)<(n); (i)++)
#define rep1(i,n) for(ll (i)=1; (i)<=(n); (i)++)
#define rep3(i,l,r) for(auto (i)=(l); (i)<(r); (i)++)

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

// �N���b�s���O  // clamp (C++17) �Ɠ���
template<class T> inline T clipped(const T& v, const T& low, const T& high) {
	return min(max(v, low), high);
}

// 2 �����x�N�g��
template<typename T> struct Vec2 {
	/*
	y �����͉�����
	x �����͉E����
	��]�͎��v��肪���iy ��������ƍl����Ɣ����v���ɂȂ�j
	*/
	T y, x;
	constexpr inline Vec2() = default;
	constexpr inline Vec2(const T& arg_y, const T& arg_x) : y(arg_y), x(arg_x) {}
	inline Vec2(const Vec2&) = default;  // �R�s�[
	inline Vec2(Vec2&&) = default;  // ���[�u
	inline Vec2& operator=(const Vec2&) = default;  // ���
	inline Vec2& operator=(Vec2&&) = default;  // ���[�u���
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
	inline Vec2 operator*(const Vec2& rhs) const {  // x + yj �Ƃ݂Ȃ�
		return Vec2(x * rhs.y + y * rhs.x, x * rhs.x - y * rhs.y);
	}
	template<typename S> inline Vec2 operator/(const S& rhs) const {
		ASSERT(rhs != 0.0, "Zero division!");
		return Vec2(y / rhs, x / rhs);
	}
	inline Vec2 operator/(const Vec2& rhs) const {  // x + yj �Ƃ݂Ȃ�
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
	inline Vec2<double> inv() const {  // x + yj �Ƃ݂Ȃ�
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
	inline double phase() const {  // [-PI, PI) �̂͂�
		return atan2(y, x);
	}
	inline double phase_deg() const {  // [-180, 180) �̂͂�
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

// ����
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

// 2 �����z��
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

// �L���[
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

// �X�^�b�N  // �R���X�g���N�^�Ăԃ^�C�~���O�Ƃ����l�����ĂȂ��ėǂ��Ȃ�
template<class T, int max_size> struct Stack {
	array<T, max_size> data;
	int right;
	inline Stack() : data(), right(0) {}
	inline Stack(const int n) : data(), right(0) { resize(n); }
	inline Stack(const int n, const T& val) : data(), right(0) { resize(n, val); }
	inline Stack(const initializer_list<T>& init) : data(), right(init.size()) {
		memcpy(&data[0], init.begin(), sizeof(T) * init.size());
	}  // ������v���H
	inline Stack(const Stack& rhs) : data(), right(rhs.right) {  // �R�s�[
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

// �n�b�V���e�[�u��  // ���܂������ł���
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
		// ���ɒl���i�[����Ă��邱�Ƃ�����
		if (key == EMPTY || key == DELETED) return Get(key + (KeyType)2);
		auto address = key & mask;
		while (data[address].first != key) address = (address + 1) & mask;
		return data[address].second;
	}
	const T& Get(const KeyType& key, const T& default_value) const {
		// �܂��l���i�[����Ă��Ȃ��ꍇ�̒l���w��
		if (key == EMPTY || key == DELETED) return Get(key + (KeyType)2, default_value);
		auto address = key & mask;
		while (true) {
			if (data[address].first == key) return data[address].second;
			else if (data[address].first == EMPTY) return default_value;
			address = (address + 1) & mask;
		}
	}
	void Set(const KeyType& key, const T& value) {
		// �܂��l���i�[����Ă��Ȃ����Ƃ�����
		if (key == EMPTY || key == DELETED) return Set(key + (KeyType)2, value);
		auto address = key & mask;
		while (data[address].first != EMPTY || data[address].first != DELETED) address = (address + 1) & mask;
		data[address].first = key;
		data[address].second = value;
	}
	T& operator[](const KeyType& key) {
		// ���݂���΂��̒l��Ԃ�
		// ���݂��Ȃ���ΐV��������ĕԂ�
		if (key == EMPTY || key == DELETED) return operator[](key + (KeyType)2);
		auto address = key & mask;
		while (data[address].first != EMPTY) {
			if (data[address].first == key) return data[address].second;
			address = (address + 1) & mask;
		}
		address = key & mask;
		while (data[address].first != EMPTY || data[address].first != DELETED) address = (address + 1) & mask;
		data[address].first = key;
		new(&data[address].second) T();
		return data[address].second;
	}
	void erase(const KeyType& key) {
		// ���ɒl���i�[����Ă��邱�Ƃ�����
		if (key == EMPTY || key == DELETED) return erase(key + (KeyType)2);
		auto address = key & mask;
		while (data[address].first != key) address = (address + 1) & mask;
		data[address].first = DELETED;
		data[address].second.~T();  // ����ǂ�����΂�����
		memset(&data[address].second, 0, sizeof(T));
	}
	void clear() {
		// �R�X�g���ł����A�f�X�g���N�^�Ƃ��Ă΂�Ȃ�
		memset(&data[0], 0, sizeof(data));
	}

};
auto a = HashMap<double>();

// ���� (�b)
inline double Time() {
	return static_cast<double>(chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now().time_since_epoch()).count()) * 1e-9;
}


// �d������
template<class VectorLike> inline void Deduplicate(VectorLike& vec) {
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
}


// 2 ���@
template<class VectorLike, typename T> inline int SearchSorted(const VectorLike& vec, const T& a) {
	return lower_bound(vec.begin(), vec.end(), a) - vec.begin();
}


// argsort
template<typename T, int n, typename result_type> inline auto Argsort(const array<T, n>& vec) {
	array<result_type, n> res;
	iota(res.begin(), res.end(), 0);
	sort(res.begin(), res.end(), [&](const result_type& l, const result_type& r) {
		return vec[l] < vec[r];
	});
	return res;
}


// popcount  // SSE 4.2 ���g���ׂ�
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

// x >> n & 1 �� 1 �ɂȂ�ŏ��� n ( x==0 �͖���` )
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

#endif  // NAGISS_LIBRARY_HPP

using ull = unsigned long long;
using i8 = int8_t;
using u8 = uint8_t;
using u16 = uint16_t;

// ����
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
			_mm256_permute2x128_si256(data, data, 0b10000001), data, 2  // _mm256_permute2x128_si256(data, data, 0b10000001) �� data �̏�ʃr�b�g���擾
		) };  // alignr(s1, s2, mask) := ((s1 �̏�ʂ� s2 �̏��) >> mask) << 16 | ((s1 �̉��ʂ� s2 �̉���) >> mask) �V�t�g�� 8 bit �P��
	}
	inline BitBoard Down() const {
		return BitBoard{ _mm256_alignr_epi8(
			data, _mm256_permute2x128_si256(data, data, 0b00001000), 16 - 2  // _mm256_permute2x128_si256(data, data, 0b00001000) �� data �̉��ʃr�b�g����ʂɎ����Ă���
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
	inline BitBoard& Expand() {  // �㉺���E�ɍL����
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
constexpr auto K = 0.02;                            // �傫���قǖ����̉��l���������Ȃ� log2/100 = 0.007 ���炢�̂Ƃ���؂̃C���t���ƒނ荇���H
const auto EXP_NEG_K = exp(-K);
auto EXP_NEG_KT = array<double, 1000>();
auto v_modified = array<double, M>();               // �^�[���ŕ␳������؂̉��l
auto NEIGHBOR = array<array<BitBoard, 16>, 256>();  // neighbor[idx][d] := idx ���狗�� d �ȓ��̏ꏊ����
auto s_begins = array<short, T + 1>();              // t ���ڂ̖�؂̍ŏ��̃C���f�b�N�X
auto order_e = array<short, M>();                   // argsort(E)
auto e_begins = array<short, T + 1>();              // order_e �̃C���f�b�N�X�ŁAt ���ڂɏ��ł����؂̍ŏ��̃C���f�b�N�X
auto start_bitboards = array<BitBoard, T>();        // ���̃^�[���ɏo�������؂̈ʒu
auto end_bitboards = array<BitBoard, T>();          // ���̃^�[���ɏ��ł����؂̈ʒu


// �r�[���T�[�`���ɕϓ�
auto t = 0;
auto future_value_table = Board<double, N, N>();  // �����������؂̉��l
auto current_value_table = Board<double, N, N>();  // �������Ă��؂̉��l
auto current_money_table = Board<short, N, N>();  // �������Ă��؂̉��l
auto current_index_table = Board<short, N, N>();  // ��؃C���f�b�N�X  // TODO: -1 �ŏ�����

void UpdateValueTable() {
	// �s������ɌĂ΂��
	
	// �o��
	rep3(idx_RCSEV, globals::s_begins[t], globals::s_begins[t + 1]) {
		const auto& rc = RC[idx_RCSEV];
		const auto& vm = v_modified[idx_RCSEV];
		const auto& v = V[idx_RCSEV];
		
		ASSERT(t == S[idx_RCSEV], "turn ������������");
		ASSERT(current_index_table.data[rc] < 0, "���ɖ�؂������");
		ASSERT(current_money_table.data[rc] == 0, "���ɖ�؂������");

		current_index_table.data[rc] = idx_RCSEV;
		current_money_table.data[rc] = v;
		current_value_table.data[rc] = vm;
		future_value_table.data[rc] -= vm;

		ASSERT(future_value_table.data[rc] >= -1e3, "�����̉��l���}�C�i�X�ɂȂ邱�Ƃ͂Ȃ��͂�����");
	}

	// ����
	rep3(idx_order, e_begins[t], e_begins[t + 1]) {
		const auto& idx_RCSEV = order_e[idx_order];

		const auto& rc = RC[idx_RCSEV];
		const auto& vm = v_modified[idx_RCSEV];
		const auto& v = V[idx_RCSEV];

		ASSERT(t == S[idx_RCSEV], "turn ������������");
		ASSERT(current_index_table.data[rc] == idx_RCSEV, "���ł������؂��Ȃ���");
		ASSERT(current_money_table.data[rc] == v, "���ł������؂��Ȃ���");

		current_index_table.data[rc] = idx_RCSEV;
		current_money_table.data[rc] = 0;
		current_value_table.data[rc] = 0.0;
	}

	t++;
}

}


struct State {
	BitBoard vegetables;
	BitBoard machines;  // ���^�[���ȍ~�͖؂�ۂ�
	short turn;  // �����u���ĂȂ���Ԃ� 0
	short n_machines;
	int money;
	double score;
	double subscore2;
	double subscore3;

	unsigned hash;  // machines �݂̂ɂ���Ĉ�ӂɒ�܂�

	struct Action {
		// �V�����u���Ƃ��� before == after �ɂ���
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
		// 1. ���n�@���ړ�������
		//   - machine ��ύX����
		//   - future_value_table �ɉ����� subscore2 �������v�Z����
		//   - (��) ��؂������Ă� money �𑝂₷�K�v�͖���
		// 2. ���̓��o�������؂ɉ����� vegetables �̃r�b�g�𗧂Ă�
		// 3. machines �� vegetables �̋��ʕ���������āA��؂����n����
		//   - ���n���������� vegetables �̃r�b�g�͐܂�
		//   - money ��������
		// 4. ��؂̏o���ɉ����� future_value_table �����������A���̃}�X�� machine ������� subscore2 �����炷
		//   - ���ۂɂ͕ϓ����������ꏊ�̃��X�g�����Atable �͐G��Ȃ�
		// 5. ���̓����ł����؂ɉ����� vegetables �̃r�b�g��܂�
		// 6. (subscore3 �͈�U�ȗ�)
		
		// �Q�Ƃ���O���̕ϐ�:
		// future_value_table
		// current_money_table
		// start_bitboards
		// end_bitboards
		// �ق�


		// Step 1
		const auto& before = action.before;
		const auto& after = action.after;
		if (before == after) {
			if (machines.Get(after)) {
				// �p�X����ꍇ
				// �������Ȃ�
			}
			else {
				// �V�����u���ꍇ
				n_machines++;
				money -= (int)n_machines * (int)n_machines * (int)n_machines;
				ASSERT(money >= 0, "money < 0");

				machines.Flip(after);
				hash ^= globals::RNT[after];
				subscore2 += globals::future_value_table.data[after];
			}
		}
		else {
			// �ړ�������ꍇ
			ASSERT(machines.Get(before), "�ړ����ɋ@�B���Ȃ���");
			ASSERT(!machines.Get(after), "�ړ���ɋ@�B�������");
			machines.Flip(before);
			machines.Flip(after);
			hash ^= globals::RNT[before] ^ globals::RNT[after];
			subscore2 += globals::future_value_table.data[after] - globals::future_value_table.data[before];
		}

		// Step 2
		vegetables |= globals::start_bitboards[turn];

		// Step 3
		auto intersection = machines & vegetables;
		vegetables ^= intersection;
		for (const auto& idx : intersection.NonzeroIndices()) {
			ASSERT(globals::current_money_table.data[idx] >= 1, "�������n���悤�Ƃ��Ă��");
			// ��ɘA�����Ă��邱�Ƃ�����
			money += n_machines * globals::current_money_table.data[idx];
		}

		// Step 4
		rep3(idx_vegetables, globals::s_begins[turn], globals::s_begins[turn + 1]) {
			const auto& idx = RC[idx_vegetables];
			const auto& vm = globals::v_modified[idx_vegetables];
			subscore2 -= machines.Get(idx) * vm;
			ASSERT(subscore2 >= -1e3, "subscore2 < 0");
		}

		// Step 5
		vegetables.data = _mm256_andnot_si256(globals::end_bitboards[turn].data, vegetables.data);


		turn++;
		if (turn == T) {
			score = money;
		}
		else {
			score = money
				+ (subscore2 + subscore3) / globals::EXP_NEG_KT[turn] * n_machines
				+ (int)((n_machines * (n_machines + 1)) / 2) * (int)((n_machines * (n_machines + 1)) / 2);  // 3 ��a
		}

		// TODO
	}
	template<class Vector>
	inline void GetNextStates(Vector& res) const {

		// TODO: n_machines �����Ȃ��Ƃ��̏���

		if (((int)n_machines + 1) * ((int)n_machines + 1) * ((int)n_machines + 1) < money) {
			// ����������Ȃ��Ȃ�ꍇ (1 ��菜��)
			if (n_machines == 1) {
				// �@�B�� 1 �̂Ƃ�
				const auto p_remove = machines.NonzeroIndices()[0];
				rep(p_add, 256) {
					if (p_remove == p_add) continue;
					auto new_state = *this;
					new_state.Do(Action{ p_remove, (u8)p_add });
					res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, (u8)p_add} });
				}
			}
			else {
				// �@�B�� 2 �ȏ�̂Ƃ�
				for (const auto& p_remove : machines.NonzeroIndices()) {
					// �t����Ȃ��������΂�
					using namespace board_index_functions;
					int neighbor_cnt = 0;  // �אڂ��鐔
					for (const auto& drul : { DownOf(p_remove), RightOf(p_remove), UpOf(p_remove), LeftOf(p_remove) }) {
						if (drul != p_remove) {
							neighbor_cnt += machines.Get(drul);
						}
					}
					if (neighbor_cnt >= 2) continue;  // p �͗t�ł͂Ȃ�

					// ���ۂɌ��炷
					auto machines_removed = machines;
					machines_removed.Flip(p_remove);

					// 1 �������@��T��
					// ���X�̘A�������� 1 �ӏ��ŗא� <=> xor �� 1 ����, 2 �y�A�� or �� xor �� 1
					auto&& cand = (machines_removed.Down() ^ machines_removed.Right() ^ machines_removed.Up() ^ machines_removed.Left())
						& ((machines_removed.Down() | machines_removed.Right()) ^ (machines_removed.Up() | machines_removed.Left()));
					for (const auto& p_add : cand.NonzeroIndices()) {
						if (p_remove == p_add) continue;
						auto new_state = *this;
						new_state.Do(Action{ p_remove, p_add });
						res.push(NewStateInfo{ new_state.score, new_state.hash, {p_remove, p_add} });
					}
				}
			}
		}
		else {
			// ����������Ă�ꍇ
			if (n_machines == 0) {
				rep(p_add, 256) {
					auto new_state = *this;
					new_state.Do(Action{ (u8)p_add, (u8)p_add });
					res.push(NewStateInfo{ new_state.score, new_state.hash, {(u8)p_add, (u8)p_add} });
				}
			}
			else {
				// 1 �������@��T��
				auto machines_copy = machines;
				auto&& cand = (machines_copy.Down() ^ machines_copy.Right() ^ machines_copy.Up() ^ machines_copy.Left())
					& ((machines_copy.Down() | machines_copy.Right()) ^ (machines_copy.Up() | machines_copy.Left()));
				for (const auto& p_add : cand.NonzeroIndices()) {
					auto new_state = *this;
					new_state.Do(Action{ p_add, p_add });
					res.push(NewStateInfo{ new_state.score, new_state.hash, {p_add, p_add} });
				}
			}
		}

		// �p�� v �� w �ɂȂ��Ƃ����l�����ׂ��H

	}

};

void Solve() {
	// ���͂��󂯎��
	{
		int buf;
		cin >> buf >> buf >> buf;
		rep(i, M) {
			cin >> R[i] >> C[i] >> S[i] >> E[i] >> V[i];
			RC[i] = ((u8)R[i] * (u8)N + (u8)C[i]);
		}
	}

	// �F�X������
	{
		using namespace globals;
		for (auto&& r : RNT) {
			r = (unsigned)rng.next();
		}
		EXP_NEG_KT[0] = 1.0;
		rep1(i, T - 1) EXP_NEG_KT[i] = EXP_NEG_KT[i - 1] * EXP_NEG_K;
		rep(i, M) {
			v_modified[i] = V[i] * EXP_NEG_KT[S[i]];
			future_value_table[{ R[i], C[i] }] += v_modified[i];
			start_bitboards[S[i]].Flip(RC[i]);
			end_bitboards[E[i]].Flip(RC[i]);
		}
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
			while (idx < S.size() && S[idx] <= turn) {
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
			while (idx < E.size() && E[order_e[idx]] <= turn) {
				idx++;
			}
			e_begins[turn + 1] = idx;
		}

		current_index_table.Fill(-1);
	}

	//globals::future_value_table.Print();
	//cout << "NEIGHBOR[200][7]" << endl;
	//globals::NEIGHBOR[200][7].Print();

	// �r�[���T�[�`
	
	{
		struct Node {
			double score;
			Node* parent_node;
			State* state;
			typename State::Action action;

			inline bool operator<(const Node& rhs) const { return score < rhs.score; }
			inline bool operator>(const Node& rhs) const { return score > rhs.score; }
		};
		static Stack<State, (int)2e6> state_buffer;
		static Stack<Node, (int)2e6> node_buffer;
		state_buffer.emplace();
		node_buffer.push({ state_buffer[0].score, nullptr, &state_buffer[0] });  // TODO �������
		static Stack<Node, 200000> q;
		Node* parent_nodes_begin = node_buffer.begin();
		Node* parent_nodes_end = node_buffer.end();
		constexpr auto hash_table_size = 14;
		constexpr int beam_width = 10;
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
						// �܂����̃n�b�V���̏�Ԃ�����
						q.push({ r.score, parent_node, nullptr, r.action });
						old_node = &q.back();
					} else if (old_node->score < r.score) {
						// �����n�b�V���ŃX�R�A�����ǂ�
						old_node->score = r.score;
						old_node->parent_node = parent_node;
						old_node->action = r.action;
					}
				}

			}
			cerr << "q.size()=" << q.size() << endl;

			if (beam_width < q.size()) {
				nth_element(q.begin(), q.begin() + beam_width, q.end(), greater<>());  // �����Ńn�b�V���e�[�u���j�󂳂�Ă���
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

		// ���ʂ��o��
		ASSERT(best_node != nullptr, "best_node ���Ȃ���");
		auto path = array<State::Action, T>();
		auto node = best_node;
		rep(i, T) {
			path[T - 1 - i] = node->action;
			node = node->parent_node;
		}
		ASSERT(node == nullptr, "���m�[�h����Ȃ���");
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


int main() {
	//test::TestBitBoard();
	Solve();

}


