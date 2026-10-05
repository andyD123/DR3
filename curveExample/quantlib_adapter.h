#pragma once
// QuantLib owns curve evaluation, dates, coupon amounts, ex-coupon rules and
// settlement conventions. This adapter snapshots their fixed-income outputs.
#include "dr3_pricing.h"
#include <ql/cashflow.hpp>
#include <ql/instruments/bonds/fixedratebond.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/settings.hpp>

namespace dr3_curve::quantlib {
using CurveHandle = QuantLib::Handle<QuantLib::YieldTermStructure>;
using BondPtr = QuantLib::ext::shared_ptr<QuantLib::FixedRateBond>;
using Serial = QuantLib::Date::serial_type;
using Table = DiscountSnapshot<Serial>;

struct Snapshot {
    std::shared_ptr<const Table> discounts;
    QuantLib::Date referenceDate;
    QuantLib::Date evaluationDate;
};

// Not an Observer: an existing snapshot intentionally remains unchanged when a
// quote or handle changes. Rebuild AND re-adapt after market/date/instrument updates.
inline Snapshot snapshotCurves(const std::vector<CurveHandle>& curves,
                               const std::vector<BondPtr>& bonds) {
    if (curves.empty() || curves.front().empty()) throw std::invalid_argument("No scenario curves");
    const auto ref = curves.front()->referenceDate();
    if (ref != QuantLib::Settings::instance().evaluationDate())
        throw std::invalid_argument("This adapter requires reference date = evaluation date");
    for (const auto& curve : curves)
        if (curve.empty() || curve->referenceDate() != ref)
            throw std::invalid_argument("Scenario curves need the same reference date");
    std::vector<Serial> dates;
    dates.push_back(ref.serialNumber());
    for (const auto& bond : bonds) {
        if (!bond) throw std::invalid_argument("Null bond");
        const auto settlement = bond->settlementDate();
        if (settlement < ref) throw std::invalid_argument("Settlement precedes curve reference date");
        dates.push_back(settlement.serialNumber());
        for (const auto& cf : bond->cashflows())
            if (cf->date() >= ref) dates.push_back(cf->date().serialNumber());
    }
    return {Table::build(std::move(dates), curves.size(), [&](Serial date, std::size_t s) {
                // Retain each curve's own interpolation and extrapolation policy.
                return curves[s]->discount(QuantLib::Date(date));
            }), ref, QuantLib::Settings::instance().evaluationDate()};
}

inline std::vector<Payment> prepareCashflows(const QuantLib::Leg& leg,
                                           QuantLib::Date cutoff,
                                           const Table& table) {
    std::vector<Payment> result;
    for (const auto& cf : leg) {
        // Explicitly exclude flows on the cutoff date, just like the benchmark
        // DiscountingBondEngine(curve, false). QuantLib supplies ex-coupon rules.
        if (!cf->hasOccurred(cutoff, false) && !cf->tradingExCoupon(cutoff))
            result.push_back({table.index(cf->date().serialNumber()), cf->amount()});
    }
    return result;
}

inline int checkedWidth(std::size_t n) {
    if (!n || n > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Invalid scenario count");
    return static_cast<int>(n);
}
inline const Table& checkedTable(const Snapshot& snapshot) {
    if (!snapshot.discounts) throw std::invalid_argument("Null snapshot");
    return *snapshot.discounts;
}

struct Prices {
    explicit Prices(std::size_t n) : npv(0.0, checkedWidth(n)),
                                    dirty(0.0, checkedWidth(n)),
                                    clean(0.0, checkedWidth(n)) {}
    Vector npv, dirty, clean;
};

class AdaptedBond {
public:
    AdaptedBond(const QuantLib::FixedRateBond& bond, const Snapshot& snapshot)
        : snapshot_(snapshot),
          npv_(snapshot.discounts, prepareCashflows(bond.cashflows(), snapshot.referenceDate, checkedTable(snapshot))),
          settlement_(snapshot.discounts, prepareCashflows(bond.cashflows(), bond.settlementDate(), checkedTable(snapshot))),
          settlementRow_(checkedTable(snapshot).index(bond.settlementDate().serialNumber())),
          notional_(bond.notional(bond.settlementDate())),
          accrued_(bond.accruedAmount(bond.settlementDate())) {
        if (QuantLib::Settings::instance().evaluationDate() != snapshot_.evaluationDate)
            throw std::invalid_argument("Evaluation date changed: rebuild snapshot");
        if (!std::isfinite(notional_) || notional_ < 0 || !std::isfinite(accrued_))
            throw std::invalid_argument("Unsupported notional/accrual");
    }
    void priceInto(Prices& out) const {
        if (out.clean.isScalar() || static_cast<std::size_t>(out.clean.size()) != snapshot_.discounts->scenarios())
            throw std::invalid_argument("Output scenario width differs");
        // Frozen valuation: no observer registrations and no QuantLib calls here.
        npv_.priceInto(out.npv);
        settlement_.priceInto(out.dirty);
        if (notional_ == 0) {
            std::fill(out.dirty.begin(), out.dirty.end(), 0.0);
            std::fill(out.clean.begin(), out.clean.end(), 0.0);
            return;
        }
        const Vector::INS scale(100.0 / notional_);
        auto atSettlement = [scale](auto value, auto df) { return (value / df) * scale; };
        ::transformM(atSettlement, out.dirty, snapshot_.discounts->row(settlementRow_));
        const Vector::INS accrued(accrued_);
        auto subtractAccrued = [accrued](auto dirty) { return dirty - accrued; };
        ::transform(subtractAccrued, out.dirty, out.clean);
    }
    Prices prices() const { Prices out(snapshot_.discounts->scenarios()); priceInto(out); return out; }
    void scalarScenario(std::size_t s, double& npv, double& dirty, double& clean) const {
        npv = npv_.scalarScenario(s);
        dirty = notional_ == 0 ? 0 : settlement_.scalarScenario(s) /
            snapshot_.discounts->row(settlementRow_)[s] * (100.0 / notional_);
        clean = notional_ == 0 ? 0 : dirty - accrued_;
    }
private:
    Snapshot snapshot_;
    PreparedLeg<Serial> npv_, settlement_;
    std::size_t settlementRow_;
    double notional_, accrued_;
};
inline AdaptedBond adaptBond(const QuantLib::FixedRateBond& bond, const Snapshot& snapshot) {
    if (!snapshot.discounts) throw std::invalid_argument("Null snapshot");
    return AdaptedBond(bond, snapshot);
}
}
