#pragma once
// Scalar dates, scalar or DR3 vector ordinates. Adapted from the DR3 Curve branch.
// Zero-rate curves use annualised continuously compounded rates and YEAR FRACTIONS.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <list>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace curve_detail {
template<class V> auto exp_value(const V& x) { using std::exp; return exp(x); }
template<class T> void finite_date(T t) {
    static_assert(std::is_arithmetic<T>::value, "Curve dates must be arithmetic");
    if (!std::isfinite(static_cast<double>(t)))
        throw std::invalid_argument("Curve date is not finite");
}
// The caller guarantees lo < t < hi. Integer differences must be formed
// before floating conversion, without signed overflow at either end of the range.
template<class T> double interpolation_weight(T t, T lo, T hi) {
    if constexpr (std::is_integral<T>::value && !std::is_same<T, bool>::value) {
        using U = std::make_unsigned_t<T>;
        const U numerator = static_cast<U>(static_cast<U>(t) - static_cast<U>(lo));
        const U denominator = static_cast<U>(static_cast<U>(hi) - static_cast<U>(lo));
        return static_cast<double>(numerator) / static_cast<double>(denominator);
    } else {
        const double left = static_cast<double>(lo), right = static_cast<double>(hi);
        const double query = static_cast<double>(t);
        const double denominator = right - left;
        if (std::isfinite(denominator)) return (query - left) / denominator;
        // Opposite-sign finite endpoints can overflow their difference. Scaling
        // each operand first keeps the fraction finite (also on MSVC).
        return (query * 0.5 - left * 0.5) / (right * 0.5 - left * 0.5);
    }
}
template<class V> std::size_t width(const V& v) {
    if constexpr (std::is_arithmetic<V>::value) {
        if (!std::isfinite(static_cast<double>(v)))
            throw std::invalid_argument("Curve value is not finite");
        return 1;
    } else {
        if (v.isScalar() || v.size() <= 0)
            throw std::invalid_argument("Curve vectors must be nonempty, nonscalar vectors");
        for (int i = 0; i < v.size(); ++i)
            if (!std::isfinite(static_cast<double>(v[i])))
                throw std::invalid_argument("Curve vector value is not finite");
        return static_cast<std::size_t>(v.size());
    }
}
template<class T, class V> void validate(const std::vector<T>& xs, const std::vector<V>& ys) {
    if (xs.empty() || xs.size() != ys.size())
        throw std::invalid_argument("Curve needs equally sized, nonempty dates and values");
    const auto n = width(ys.front());
    for (std::size_t i = 0; i < xs.size(); ++i) {
        finite_date(xs[i]);
        if (i && !(xs[i - 1] < xs[i]))
            throw std::invalid_argument("Curve dates must be strictly increasing");
        if (width(ys[i]) != n) throw std::invalid_argument("Curve scenario widths differ");
    }
}
template<class T, class V> void check_query(T t, const std::vector<T>& xs, const std::vector<V>& ys) {
    finite_date(t);
    if (xs.empty() || xs.size() != ys.size()) throw std::logic_error("Curve is not initialised");
}
}

template<class T, class V> struct Extrap_end_flat_last {
    static V extrap(T, const std::vector<V>& ys, const std::vector<T>& xs) {
        if (xs.empty() || ys.empty()) throw std::logic_error("Empty curve");
        return ys.back();
    }
};
template<class T, class V> struct Extrap_start_flat_first {
    static V extrap(T, const std::vector<V>& ys, const std::vector<T>& xs) {
        if (xs.empty() || ys.empty()) throw std::logic_error("Empty curve");
        return ys.front();
    }
};

template<class T, class V> struct LinearInterp {
    static V calc(T t, const std::vector<V>& ys, const std::vector<T>& xs) {
        curve_detail::check_query(t, xs, ys);
        const auto hi = std::lower_bound(xs.begin(), xs.end(), t);
        if (hi != xs.end() && *hi == t) return ys[static_cast<std::size_t>(hi - xs.begin())];
        if (hi == xs.begin() || hi == xs.end()) throw std::out_of_range("Interpolation outside pillars");
        const auto i = static_cast<std::size_t>(hi - xs.begin());
        const double w = curve_detail::interpolation_weight(t, xs[i - 1], xs[i]);
        return ys[i - 1] + (ys[i] - ys[i - 1]) * w;
    }
    static V start(T t, const std::vector<V>& ys, const std::vector<T>& xs) {
        return Extrap_start_flat_first<T,V>::extrap(t, ys, xs);
    }
    static V end(T t, const std::vector<V>& ys, const std::vector<T>& xs) {
        return Extrap_end_flat_last<T,V>::extrap(t, ys, xs);
    }
};
template<class T, class V> struct FlatInterp : LinearInterp<T,V> {
    static V calc(T t, const std::vector<V>& ys, const std::vector<T>& xs) {
        curve_detail::check_query(t, xs, ys);
        const auto hi = std::upper_bound(xs.begin(), xs.end(), t);
        if (hi == xs.begin() || xs.back() < t) throw std::out_of_range("Interpolation outside pillars");
        return ys[static_cast<std::size_t>(hi - xs.begin() - 1)];
    }
};
template<class T, class V> struct ZeroInterp {
    static V discount(T t, const V& r) {
        curve_detail::finite_date(t);
        if (t < T(0)) throw std::out_of_range("Discount-curve time is negative");
        return curve_detail::exp_value(r * (-static_cast<double>(t)));
    }
    static V calc(T t, const std::vector<V>& ys, const std::vector<T>& xs) {
        return discount(t, LinearInterp<T,V>::calc(t, ys, xs));
    }
    // Extrapolate ZERO RATES, not raw rates mistaken for discount factors.
    static V start(T t, const std::vector<V>& ys, const std::vector<T>&) { return discount(t, ys.front()); }
    static V end(T t, const std::vector<V>& ys, const std::vector<T>&) { return discount(t, ys.back()); }
};

template<class T, class V, class INTERP_POLICY = LinearInterp<T,V>,
         class EXTRAP_END_POLICY = void, class EXTRAP_START_POLICY = void>
class Curve {
public:
    template<class IT_T, class IT_V> void setValues(IT_T xb, IT_T xe, IT_V yb, IT_V ye) {
        std::vector<T> xs(xb, xe);
        std::vector<V> ys(yb, ye);
        curve_detail::validate(xs, ys);
        dates_.swap(xs); values_.swap(ys); // Failed validation leaves the old curve intact.
    }
    V valueAt(T t) const {
        curve_detail::check_query(t, dates_, values_);
        if (t < dates_.front()) {
            if constexpr (std::is_void<EXTRAP_START_POLICY>::value) return INTERP_POLICY::start(t, values_, dates_);
            else return EXTRAP_START_POLICY::extrap(t, values_, dates_);
        }
        if (dates_.back() < t) {
            if constexpr (std::is_void<EXTRAP_END_POLICY>::value) return INTERP_POLICY::end(t, values_, dates_);
            else return EXTRAP_END_POLICY::extrap(t, values_, dates_);
        }
        return INTERP_POLICY::calc(t, values_, dates_);
    }
private:
    std::vector<T> dates_;
    std::vector<V> values_;
};

// References survive hits and unrelated insertions, but NOT eviction, replacement,
// clear/reset, assignment or destruction. A mutable cache is not thread safe.
template<class K, class V> class lru_cache {
    using Item = std::pair<K,V>;
    using List = std::list<Item>;
    using Iterator = typename List::iterator;
public:
    explicit lru_cache(std::size_t capacity) : capacity_(capacity) {
        if (!capacity) throw std::invalid_argument("LRU capacity must be positive");
    }
    lru_cache(const lru_cache& other) : capacity_(other.capacity_), items_(other.items_) { reindex(); }
    lru_cache& operator=(const lru_cache& other) {
        if (this != &other) { lru_cache copy(other); swap(copy); }
        return *this;
    }
    lru_cache(lru_cache&&) = default;
    lru_cache& operator=(lru_cache&&) = default;
    void swap(lru_cache& other) {
        std::swap(capacity_, other.capacity_); items_.swap(other.items_); index_.swap(other.index_);
    }
    const V* find(const K& key) {
        const auto it = index_.find(key);
        if (it == index_.end()) return nullptr;
        items_.splice(items_.begin(), items_, it->second);
        return &it->second->second;
    }
    const V& get(const K& key) {
        const auto* p = find(key);
        if (!p) throw std::out_of_range("Key is not cached");
        return *p;
    }
    const V& put(const K& key, V value) {
        const auto it = index_.find(key);
        if (it != index_.end()) {
            it->second->second = std::move(value);
            items_.splice(items_.begin(), items_, it->second);
            return it->second->second;
        }
        items_.emplace_front(key, std::move(value));
        try { index_.emplace(items_.front().first, items_.begin()); }
        catch (...) { items_.pop_front(); throw; }
        if (items_.size() > capacity_) { index_.erase(items_.back().first); items_.pop_back(); }
        return items_.front().second;
    }
    bool exists(const K& key) const { return index_.find(key) != index_.end(); }
    std::size_t size() const { return items_.size(); }
    void clear() { index_.clear(); items_.clear(); }
private:
    void reindex() { for (auto it = items_.begin(); it != items_.end(); ++it) index_.emplace(it->first, it); }
    std::size_t capacity_;
    List items_;
    std::unordered_map<K,Iterator> index_;
};

// Compatibility spelling. Cached and uncached zero interpolation now have the
// SAME mathematical definition; the old undocumented MyCalc recurrence is gone.
template<class T, class V> struct ZeroInterpCached : ZeroInterp<T,V> {};

template<class T, class V, class INTERP_POLICY = LinearInterp<T,V>,
         class EXTRAP_END_POLICY = void, class EXTRAP_START_POLICY = void>
class Curve2 {
public:
    explicit Curve2(std::size_t capacity) : cache_(capacity) {}
    template<class IT_T, class IT_V> void setValues(IT_T xb, IT_T xe, IT_V yb, IT_V ye) {
        curve_.setValues(xb, xe, yb, ye);
        cache_.clear(); // Invalidate interpolation AND extrapolation results.
    }
    const V& valueAtRef(T t) {
        curve_detail::finite_date(t);
        if (const auto* p = cache_.find(t)) return *p;
        return cache_.put(t, curve_.valueAt(t));
    }
    V valueAt(T t) { return valueAtRef(t); } // Original owning/copying interface.
    void clearCache() { cache_.clear(); }
    std::size_t cachedDates() const { return cache_.size(); }
private:
    Curve<T,V,INTERP_POLICY,EXTRAP_END_POLICY,EXTRAP_START_POLICY> curve_;
    lru_cache<T,V> cache_;
};
