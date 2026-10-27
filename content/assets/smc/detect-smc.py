"""Chronological SMC labels on completed OHLC bars; not a trading simulator.

Usage: python3 detect-smc.py example.json
Only the Python standard library is required. Bar numbers are one-based.
Strict pivots use k bars on either side. Sweep/break levels must have been
confirmed BEFORE the current bar; pivots confirmed now apply to later bars.
"""

import json
import math
from pathlib import Path
import sys


def detect(bars, k=2):
    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        raise ValueError("k must be a positive integer")
    events = []
    known_low = known_high = None
    broken_highs = set()

    def emit(kind, t, **values):
        events.append({"kind": kind, "known_at": t + 1, **values})

    for t, bar in enumerate(bars):
        o, h, l, c = (bar[key] for key in ("open", "high", "low", "close"))
        if not all(math.isfinite(v) for v in (o, h, l, c)):
            raise ValueError(f"Nonfinite price at bar {t+1}")
        if not (l <= min(o, c) <= max(o, c) <= h):
            raise ValueError(f"Invalid OHLC at bar {t+1}")

        # Levels were known before this bar, not confirmed using this bar.
        if known_low and l < known_low[1] < c:
            emit("low_sweep", t, reference_bar=known_low[0],
                 level=known_low[1], low=l, close=c)
        if known_high and known_high[0] not in broken_highs and c > known_high[1]:
            emit("upward_break", t, reference_bar=known_high[0],
                 level=known_high[1], close=c)
            broken_highs.add(known_high[0])
            for j in range(t-1, max(-1, t-6), -1):
                if bars[j]["close"] < bars[j]["open"]:
                    emit("ob_candidate", t, origin_bar=j+1,
                         zone=[bars[j]["low"], bars[j]["high"]])
                    break

        if t >= 2:
            if l > bars[t-2]["high"]:
                emit("bullish_fvg", t, first_bar=t-1,
                     zone=[bars[t-2]["high"], l])
            if h < bars[t-2]["low"]:
                emit("bearish_fvg", t, first_bar=t-1,
                     zone=[h, bars[t-2]["low"]])

        i = t - k
        if i >= k:
            neighbours = bars[i-k:i] + bars[i+1:t+1]
            if all(bars[i]["low"] < b["low"] for b in neighbours):
                known_low = (i+1, bars[i]["low"])
                emit("swing_low", t, origin_bar=i+1, price=known_low[1])
            if all(bars[i]["high"] > b["high"] for b in neighbours):
                known_high = (i+1, bars[i]["high"])
                emit("swing_high", t, origin_bar=i+1, price=known_high[1])
    return events


def check_example(document):
    bars = document["bars"]
    events = detect(bars)
    # Each truncated run must equal the already available full-run events.
    for n in range(len(bars) + 1):
        assert detect(bars[:n]) == [e for e in events if e["known_at"] <= n]
    assert any(e == {"kind": "swing_low", "known_at": 5,
                     "origin_bar": 3, "price": 98} for e in events)
    assert any(e == {"kind": "swing_high", "known_at": 8,
                     "origin_bar": 6, "price": 107} for e in events)
    assert any(e["kind"] == "low_sweep" and e["known_at"] == 9
               and e["level"] == 98 for e in events)
    assert any(e["kind"] == "upward_break" and e["known_at"] == 11
               and e["level"] == 107 for e in events)
    assert any(e["kind"] == "ob_candidate" and e["known_at"] == 11
               and e["origin_bar"] == 8 and e["zone"] == [99.5, 104]
               for e in events)
    assert [e for e in events if e["kind"] == "bullish_fvg"] == [
        {"kind": "bullish_fvg", "known_at": 12, "first_bar": 10,
         "zone": [102, 105]}]
    assert any(e["kind"] == "bearish_fvg" and e["known_at"] == 9
               and e["zone"] == [101, 102] for e in events)

    # Equality is not a strict pivot or a gap; reject invalid bars.
    flat = [{"open": 1, "high": 1, "low": 1, "close": 1}] * 8
    assert detect(flat) == []
    try:
        detect([{"open": 2, "high": 1, "low": 0, "close": 1}])
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid OHLC accepted")

    # Independent check of the separate binary-payoff arithmetic.
    risk = document["separate_risk_example"]
    distance = risk["entry"] - risk["stop"]
    win = (risk["target"] - risk["entry"] - risk["all_in_cost_per_unit"]) / distance
    loss = (risk["stop"] - risk["entry"] - risk["all_in_cost_per_unit"]) / distance
    p_break_even = -loss / (win - loss)
    assert math.isclose(win, 1.95) and math.isclose(loss, -1.05)
    assert math.isclose(p_break_even, .35)
    for p, expected in ((.30, -.15), (.35, 0), (.40, .15)):
        assert math.isclose(p*win + (1-p)*loss, expected, abs_tol=1e-12)
    return events


if __name__ == "__main__":
    filename = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("example.json")
    document = json.loads(filename.read_text(encoding="utf-8"))
    events = check_example(document)
    print(json.dumps({"checks": "passed", "historical_backtest": False,
                      "events": events}, ensure_ascii=False, indent=2))
