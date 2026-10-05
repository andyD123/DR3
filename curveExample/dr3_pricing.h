#pragma once
#include "../Vectorisation/VecX/dr3.h"
#include "curve.h"
#include <limits>
#include <memory>
#include <vector>

namespace dr3_curve {
#if defined(DR3_CURVE_AVX512)
using Vector = DRC::VecD8D::VecXX;
inline constexpr const char* backend_name = "AVX512 / 8 doubles";
#elif defined(DR3_CURVE_SSE2)
using Vector = DRC::VecD2D::VecXX;
inline constexpr const char* backend_name = "SSE2 / 2 doubles";
#else
using Vector = DRC::VecD4D::VecXX;
inline constexpr const char* backend_name = "AVX2 / 4 doubles";
#endif

// One independently owned contiguous scenario row per unique date. Immutable
// snapshots avoid both LRU eviction lifetimes and invisible stale-cache updates.
// Evaluator is scalar(date, scenario): it can wrap ANY existing curve API.
template<class Key> class DiscountSnapshot {
public:
    template<class Evaluator>
    static std::shared_ptr<const DiscountSnapshot> build(std::vector<Key> dates,
                                                        std::size_t scenarios,
                                                        Evaluator evaluate) {
        if (!scenarios || scenarios > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument("Invalid scenario count");
        for (const auto& t : dates) curve_detail::finite_date(t);
        std::sort(dates.begin(), dates.end());
        dates.erase(std::unique(dates.begin(), dates.end()), dates.end());
        auto table = std::shared_ptr<DiscountSnapshot>(new DiscountSnapshot);
        table->dates_ = std::move(dates);
        table->scenarios_ = scenarios;
        table->rows_.reserve(table->dates_.size());
        for (const auto& t : table->dates_) {
            Vector row(0.0, static_cast<int>(scenarios));
            for (std::size_t s = 0; s < scenarios; ++s) {
                const double df = evaluate(t, s);
                if (!(df > 0) || !std::isfinite(df))
                    throw std::domain_error("Discount factor must be finite and positive");
                row[s] = df;
            }
            table->rows_.push_back(std::move(row));
        }
        return table;
    }
    std::size_t scenarios() const { return scenarios_; }
    std::size_t dates() const { return dates_.size(); }
    std::size_t index(Key t) const {
        curve_detail::finite_date(t);
        const auto it = std::lower_bound(dates_.begin(), dates_.end(), t);
        if (it == dates_.end() || *it != t) throw std::out_of_range("Date is not in this snapshot");
        return static_cast<std::size_t>(it - dates_.begin());
    }
    const Vector& row(std::size_t i) const { return rows_.at(i); }
private:
    DiscountSnapshot() = default;
    std::size_t scenarios_ = 0;
    std::vector<Key> dates_;
    std::vector<Vector> rows_;
};

struct Payment { std::size_t row; double amount; };

template<class Key> class PreparedLeg {
public:
    using Snapshot = DiscountSnapshot<Key>;
    PreparedLeg(std::shared_ptr<const Snapshot> snapshot, std::vector<Payment> payments)
        : snapshot_(std::move(snapshot)), payments_(std::move(payments)) {
        if (!snapshot_) throw std::invalid_argument("Null discount snapshot");
        for (const auto& cf : payments_) {
            if (cf.row >= snapshot_->dates() || !std::isfinite(cf.amount))
                throw std::invalid_argument("Invalid prepared cash flow");
        }
    }
    void priceInto(Vector& prices) const {
        if (prices.isScalar() || static_cast<std::size_t>(prices.size()) != snapshot_->scenarios())
            throw std::invalid_argument("Output scenario width differs");
        std::fill(prices.begin(), prices.end(), 0.0);
        for (const auto& cf : payments_) {
            const Vector::INS amount(cf.amount);
            auto discountAndAccumulate = [amount](auto pv, auto df) {
                return mul_add(amount, df, pv);
            };
            ::transformM(discountAndAccumulate, prices, snapshot_->row(cf.row));
        }
    }
    Vector prices() const {
        Vector out(0.0, static_cast<int>(snapshot_->scenarios()));
        priceInto(out);
        return out;
    }
    double scalarScenario(std::size_t s) const {
        if (s >= snapshot_->scenarios()) throw std::out_of_range("Scenario index");
        double pv = 0;
        for (const auto& cf : payments_) pv += cf.amount * snapshot_->row(cf.row)[s];
        return pv;
    }
private:
    std::shared_ptr<const Snapshot> snapshot_;
    std::vector<Payment> payments_;
};
}
