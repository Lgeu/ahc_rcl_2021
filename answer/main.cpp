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
//#pragma GCC target("sse4")
//#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
//#pragma GCC optimize("O3")
//#pragma GCC optimize("Ofast")
//#pragma GCC optimize("unroll-loops")
#endif

// ========================== macroes ==========================

#define rep(i,n) for(ll (i)=0; (i)<(n); (i)++)
#define rep1(i,n) for(ll (i)=1; (i)<=(n); (i)++)
#define rep3(i,s,n) for(ll (i)=(s); (i)<(n); (i)++)

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

// 盤
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
	inline Stack(initializer_list<T> init) :
		data(init.begin(), init.end()), right(init.size()) {}  // これ長さが最大じゃないと動かない
	inline Stack(const Stack& rhs) : data(), right(rhs.right) {  // コピー
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
	template <class... Args> inline void emplace(const Args&... args) {
		ASSERT_RANGE(right, 0, max_size);
		data[right] = T(args...);
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
};


// 時間 (秒)
inline double time() {
	return static_cast<double>(chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now().time_since_epoch()).count()) * 1e-9;
}


// 重複除去
template<typename T> inline void deduplicate(vector<T>& vec) {
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
}


template<typename T> inline int search_sorted(const vector<T>& vec, const T& a) {
	return lower_bound(vec.begin(), vec.end(), a) - vec.begin();
}


// popcount  // SSE 4.2 を使うべき
inline int popcount(const unsigned int& x) {
#ifdef _MSC_VER
	return (int)__popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}
inline int popcount(const unsigned long long& x) {
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

#endif  // NAGISS_LIBRARY_HPP

using ull = unsigned long long;
using i8 = int8_t;
using u8 = uint8_t;
using u16 = uint16_t;

// 入力
constexpr int N = 16;
constexpr int M = 5000;
constexpr int T = 1000;
array<int, M> R;
array<int, M> C;
array<int, M> S;
array<int, M> E;
array<int, M> V;


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
	


	void Print() const {
		for (const auto& row : Rows()) {
			rep(i, 16) cout << (row >> i & 1) << " \n"[i == 15];
		}
	}
};

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
auto rng = Random(42);
auto rnt = array<unsigned, 10000>();  // random number table
constexpr auto K = 0.01;  // 大きいほど未来の価値が小さくなる log2/100 = 0.007 くらいのとき野菜のインフレと釣り合う？
auto exp_neg_Kt = array<double, 1000>();
auto exp_neg_K = exp(-K);
auto v_modified = array<double, M>();

// ビームサーチ中に変動
auto t = 0;
auto future_value_table = Board<double, N, N>();  // 将来生える野菜の価値
auto current_value_table = Board<double, N, N>();  // 今生えてる野菜の価値

}


struct State {
	BitBoard vegetables;
	BitBoard machines;  // 一定ターン以降は木を保つ
	short turn;
	short n_machines;
	int money;
	double score;
	unsigned hash;

	struct Action {
		// 新しく置くときは before == after にする
		u8 before, after;
	};
	struct NewStateInfo {
		double score;
		Action action;
		unsigned hash;
	};
	inline void Do(const Action& action) {
		const auto& before = action.before;
		const auto& after = action.after;
		if (before == after) {
			if (machines.Get(after)) {
				// パス
				// TODO
			}
			else {
				// 新しく置く
				n_machines++;
				money -= n_machines * n_machines * n_machines;
				ASSERT(money >= 0, "money < 0");

				// TODO
			}
		}
		else {
			// 移動させる

		}

		score = 

		// TODO
	}
	template<class Vector>
	inline void GetNextStates(Vector& res) {

		auto aaa = []() {
			// 1 個足す方法を探す処理
		};

		auto put = []() {
			// 1 個足したときのスコアとハッシュを計算する処理
		};

		// 資金が足りないなら葉を 1 個取り除く
		for (const auto& p : machines.NonzeroIndices()) {
			// 葉じゃなかったら飛ばす
			using namespace board_index_functions;
			int cnt = 0;  // 隣接する数
			for (const auto& drul : { DownOf(p), RightOf(p), UpOf(p), LeftOf(p) }) {
				if (drul != p) {
					cnt += machines.Get(drul);
				}
			}
			if (cnt >= 2) continue;  // p は葉ではない

			// 実際に減らす

			// 1 個足す方法を探す

			// 実際に減らしていたところを戻す

		}
		
		// 1 個足す

		// 角で v が w になるやつ

		// TODO
	}


};


void Solve() {
	// 入力を受け取る
	{
		int buf;
		cin >> buf >> buf >> buf;
		rep(i, M) cin >> R[i] >> C[i] >> S[i] >> E[i] >> V[i];
	}

	// 色々初期化
	{
		using namespace globals;
		for (auto&& r : rnt) {
			r = (unsigned)rng.next();
		}
		exp_neg_Kt[0] = 1.0;
		rep1(i, T - 1) exp_neg_Kt[i] = exp_neg_Kt[i - 1] * exp_neg_K;
		rep(i, M) {
			v_modified[i] = V[i] * exp_neg_Kt[S[i]];
			future_value_table[{ R[i], C[i] }] += v_modified[i];
		}

	}


	// ビームサーチ
	{
		struct Candidate {
			double score;
			State* ptr_parent_state;
			typename State::Action action;
			inline Candidate() : score(0.0), ptr_parent_state(nullptr), action() {}
			inline Candidate(const double& a_score, State* const a_ptr_parent_state, const typename State::Action a_action) :
				score(a_score), ptr_parent_state(a_ptr_parent_state), action(a_action) {}
			inline bool operator<(const Candidate& rhs) const { return score < rhs.score; }
			inline bool operator>(const Candidate& rhs) const { return score > rhs.score; }
		};
		static Stack<State, (int)1e7> state_buffer;
		state_buffer.push();  // TODO 初期状態
		static Stack<Candidate, 200000> candidates;
		static Stack<typename State::NewStateInfo, 10000> next_state_buffer;
		State* parent_states_begin = state_buffer.begin();
		State* parent_states_end = state_buffer.end();
		static bitset<1 << hash_table_size> candidates_contains;  // 2^26 bits == 8 MB
		

		constexpr int beam_width = 1000;
		rep(t, T) {
			for(auto parent_state = parent_states_begin; parent_state != parent_states_end; parent_state++){
				parent_state->GetNextStates(next_state_buffer);
				for(const auto& r : next_state_buffer){
					if (!candidates_contains[r.hash & (1u << hash_table_size) - 1u]) {
						// まだそのハッシュの状態が無い
						candidates.emplace(r.score, &parent_state, r.action);
						candidates_contains[r.hash & (1u << hash_table_size) - 1u] = true;
					}
				}
				next_state_buffer.clear();
			}
			if (beam_width < candidates.size()) {

			}
		}

	}

	


}



int main() {
	Solve();

}


