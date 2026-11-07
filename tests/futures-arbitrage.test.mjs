import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { calculateCalendarSpread } from '../content/assets/futures-arbitrage/calculator.js';

const example = { nearOpen: 4000, farOpen: 4100, nearClose: 4050, farClose: 4110 };
const ledger = JSON.parse(readFileSync(new URL('../content/assets/futures-arbitrage/example.json', import.meta.url), 'utf8'));
const nearlyEqual = (actual, expected, message) => assert.ok(Math.abs(actual - expected) < 1e-9, message ?? `${actual} ≠ ${expected}`);

test('a narrowing spread earns 400 gross and 360 after both-leg costs', () => {
  assert.deepEqual(calculateCalendarSpread(example), {
    nearPnl: 500,
    farPnl: -100,
    gross: 400,
    net: 360,
    spreadOpen: 100,
    spreadClose: 60,
    breakEvenFar: 4146,
    quantity: 10,
  });
});

test('a widening spread loses money for long-near and earns money for short-near', () => {
  const widening = { ...example, farClose: 4190 };
  const longNear = calculateCalendarSpread(widening);
  const shortNear = calculateCalendarSpread({ ...widening, longNear: false });
  assert.equal(longNear.gross, -400);
  assert.equal(longNear.net, -440);
  assert.equal(shortNear.gross, 400);
  assert.equal(shortNear.net, 360);
});

test('a common closing-price move changes each leg but leaves spread P&L unchanged', () => {
  const base = calculateCalendarSpread(example);
  for (const shift of [-5000, -200, 200, 5000]) {
    const shifted = calculateCalendarSpread({
      ...example,
      nearClose: example.nearClose + shift,
      farClose: example.farClose + shift,
    });
    assert.equal(shifted.nearPnl, base.nearPnl + shift * base.quantity);
    assert.equal(shifted.farPnl, base.farPnl - shift * base.quantity);
    assert.equal(shifted.gross, base.gross);
    assert.equal(shifted.net, base.net);
  }
});

test('reversing the position reverses gross P&L but still incurs costs', () => {
  const longNear = calculateCalendarSpread(example);
  const shortNear = calculateCalendarSpread({ ...example, longNear: false });
  assert.equal(shortNear.nearPnl, -longNear.nearPnl);
  assert.equal(shortNear.farPnl, -longNear.farPnl);
  assert.equal(shortNear.gross, -longNear.gross);
  assert.equal(shortNear.net, -440);
  assert.equal(longNear.net + shortNear.net, -80);
  assert.equal(shortNear.breakEvenFar, 4154);
});

test('calculated break-even price covers costs in either direction', () => {
  for (const longNear of [true, false]) {
    for (const cost of [0, 40, 240]) {
      const inputs = { ...example, longNear, multiplier: 5, lots: 2, cost };
      const { breakEvenFar } = calculateCalendarSpread(inputs);
      const breakEven = calculateCalendarSpread({ ...inputs, farClose: breakEvenFar });
      assert.equal(breakEven.net, 0);
      assert.equal(breakEven.gross, cost);

      const onePriceUnitHigher = calculateCalendarSpread({ ...inputs, farClose: breakEvenFar + 1 });
      assert.equal(onePriceUnitHigher.net, longNear ? -10 : 10);
    }
  }
});

test('changing common contract multiplier and lots scales exposure, not the entered total cost', () => {
  const larger = calculateCalendarSpread({ ...example, multiplier: 20, lots: 3 });
  assert.equal(larger.quantity, 60);
  assert.equal(larger.nearPnl, 3000);
  assert.equal(larger.farPnl, -600);
  assert.equal(larger.gross, 2400);
  assert.equal(larger.net, 2360);
  const proportionalCost = calculateCalendarSpread({ ...example, multiplier: 20, lots: 3, cost: 240 });
  assert.equal(proportionalCost.net, 6 * 360);
});

test('zero or negative prices remain valid; unchanged spreads still incur costs', () => {
  const result = calculateCalendarSpread({ nearOpen: -30, farOpen: -20, nearClose: 0, farClose: 10 });
  assert.equal(result.nearPnl, 300);
  assert.equal(result.farPnl, -300);
  assert.equal(result.gross, 0);
  assert.equal(result.net, -40);
});

test('required prices and all numeric parameters reject nonnumeric or nonfinite input', () => {
  assert.throws(() => calculateCalendarSpread(), RangeError);
  for (const key of ['nearOpen', 'farOpen', 'nearClose', 'farClose', 'multiplier', 'lots', 'cost']) {
    for (const value of [NaN, Infinity, -Infinity, '', '10', null, true]) {
      assert.throws(() => calculateCalendarSpread({ ...example, [key]: value }), RangeError, `${key}: ${String(value)}`);
    }
  }
});

test('multiplier, lot count and total cost enforce their independent domain constraints', () => {
  for (const multiplier of [0, -1]) {
    assert.throws(() => calculateCalendarSpread({ ...example, multiplier }), RangeError);
  }
  for (const lots of [0, -1, 0.5, 1.5]) {
    assert.throws(() => calculateCalendarSpread({ ...example, lots }), RangeError);
  }
  assert.throws(() => calculateCalendarSpread({ ...example, cost: -1 }), RangeError);
  assert.equal(calculateCalendarSpread({ ...example, multiplier: 0.5, cost: 0 }).net, 20);
});

test('direction accepts only booleans instead of coercing strings or numeric flags', () => {
  for (const longNear of ['true', 'false', '', 0, 1, null]) {
    assert.throws(() => calculateCalendarSpread({ ...example, longNear }), RangeError);
  }
});

test('closing the spot and futures legs realizes the change in basis, not the initial delivery spread', () => {
  const carry = ledger.cash_and_carry;
  const exits = carry.basis_exits;
  assert.equal(exits.entry_basis_spot_minus_futures, -60);
  assert.deepEqual(exits.scenarios.map(s => s.exit_basis_spot_minus_futures), [0, -20, -50]);
  assert.deepEqual(exits.scenarios.map(s => s.gross_pnl), [600, 400, 100]);
  assert.deepEqual(exits.scenarios.map(s => s.net_pnl), [300, 100, -200]);
  for (const exit of exits.scenarios) {
    const totalSales = exit.spot_exit * carry.quantity + (carry.futures_entry - exit.futures_exit) * carry.quantity;
    const cashPaid = carry.spot_entry * carry.quantity + exit.illustrative_fixed_cost;
    assert.equal(totalSales - cashPaid, exit.net_pnl);
    assert.equal(exit.gross_pnl, carry.quantity * (exit.exit_basis_spot_minus_futures - exits.entry_basis_spot_minus_futures));
  }
  // With 30 CNY/tonne fixed costs, -30 is the zero-net exit basis.
  assert.equal((-30 - exits.entry_basis_spot_minus_futures) * carry.quantity - carry.assumed_complete_carrying_cost, 0);
});

test('the 300 CNY carry cost finances the entire allocation and reconciles to a 2030 delivery break-even', () => {
  const carry = ledger.cash_and_carry;
  const c = carry.cost_breakdown;
  assert.equal(c.total_funded_capital, 25000);
  assert.equal(c.spot_purchase_capital + c.futures_and_reserve_capital, c.total_funded_capital);
  nearlyEqual(c.spot_purchase_capital * c.annual_simple_financing_rate * c.holding_days / c.day_count_basis, c.spot_financing_cost);
  nearlyEqual(c.futures_and_reserve_capital * c.annual_simple_financing_rate * c.holding_days / c.day_count_basis, c.futures_and_reserve_financing_cost);
  assert.deepEqual([c.spot_financing_cost, c.futures_and_reserve_financing_cost, c.total_financing_cost], [120, 30, 150]);
  assert.equal(c.storage_per_tonne_per_day * carry.quantity * c.holding_days, c.storage_cost);
  assert.equal(c.total_financing_cost + c.storage_cost + c.other_complete_costs + c.additional_net_tax_cost, c.total_cost);
  assert.equal(c.total_cost, 300);
  assert.equal(c.delivery_break_even_futures_price, 2030);
  assert.equal((c.delivery_break_even_futures_price - carry.spot_entry) * carry.quantity - c.total_cost, 0);
});

test('the 540 CNY margin top-up moves existing reserve without inventing cash or charging it twice', () => {
  const carry = ledger.cash_and_carry;
  const f = carry.funding_ledger;
  const lastMark = carry.interim_path.at(-1);
  assert.equal(f.initial_spot_purchase + f.initial_futures_account_cash + f.initial_external_reserve, f.initial_borrowed_capital);
  assert.equal(f.initial_external_reserve, 1940);
  assert.equal(f.transfer_from_reserve_to_account, -lastMark.available_funds);
  assert.equal(lastMark.futures_account_equity + f.transfer_from_reserve_to_account, f.account_equity_after_top_up);
  assert.equal(f.initial_external_reserve - f.transfer_from_reserve_to_account, f.remaining_external_reserve);
  assert.equal(f.account_equity_after_top_up, lastMark.required_margin);
  assert.equal(f.account_equity_after_top_up + f.remaining_external_reserve,
    lastMark.futures_account_equity + f.initial_external_reserve);
  assert.equal(f.remaining_external_reserve, 1400);
});

test('expiry cash, futures P&L, expenses and principal independently reconcile to 300 CNY profit', () => {
  const carry = ledger.cash_and_carry;
  const f = carry.funding_ledger;
  assert.equal(f.physical_delivery_proceeds, f.assumed_delivery_price * carry.quantity);
  assert.equal(f.final_futures_cumulative_pnl, (carry.futures_entry - f.assumed_delivery_price) * carry.quantity);
  assert.equal(f.futures_account_cash_released, f.initial_futures_account_cash + f.transfer_from_reserve_to_account + f.final_futures_cumulative_pnl);
  assert.equal(f.futures_account_cash_released, 3200);
  assert.equal(f.total_cash_before_costs, f.physical_delivery_proceeds + f.futures_account_cash_released + f.remaining_external_reserve);
  assert.equal(f.total_cash_before_costs, 25600);
  assert.equal(f.total_cash_before_costs - f.total_costs_paid_at_end, f.cash_after_costs);
  assert.equal(f.cash_after_costs, 25300);
  assert.equal(f.loan_principal_repaid, f.initial_borrowed_capital);
  assert.equal(f.cash_after_costs - f.loan_principal_repaid, f.net_cash_after_repayment);
  assert.equal(f.net_cash_after_repayment, 300);
  assert.equal(f.net_cash_after_repayment, carry.assumed_delivery_net_profit);
  assert.equal(f.net_profit_over_capital_committed_for_60_days, 0.012);
  assert.equal(f.net_profit_over_capital_committed_for_60_days, f.net_cash_after_repayment / f.initial_borrowed_capital);
});

test('calendar holding costs and exit scenarios distinguish a thesis from the zero-profit condition', () => {
  const spread = ledger.calendar_spread;
  const c = spread.carry_reference;
  nearlyEqual(c.financing_per_tonne, spread.entry.near * c.annual_simple_financing_rate * c.delivery_gap_days / 365);
  assert.equal(c.storage_per_tonne, c.storage_per_tonne_per_day * c.delivery_gap_days);
  assert.equal(c.base_carry_per_tonne, c.financing_per_tonne + c.storage_per_tonne + c.other_verified_costs_per_tonne);
  assert.deepEqual([c.base_carry_per_tonne, c.stress_carry_per_tonne, c.break_even_exit_spread], [40, 60, 96]);
  assert.deepEqual(spread.exit_scenarios.map(s => s.net_pnl), [560, 360, -40, -640]);
  for (const scenario of spread.exit_scenarios) {
    const calculation = calculateCalendarSpread({
      nearOpen: spread.entry.near, farOpen: spread.entry.far,
      nearClose: 4000, farClose: 4000 + scenario.exit_far_minus_near,
      multiplier: spread.quantity_each_leg, lots: 1, cost: spread.assumed_total_trading_cost,
    });
    assert.equal(calculation.gross, scenario.gross_pnl);
    assert.equal(calculation.net, scenario.net_pnl);
  }
});

test('the statistical hedge uses contract multipliers and reproduces both successful and adverse leg cash flows', () => {
  const s = ledger.statistical_spread;
  const p = s.positions;
  const params = s.parameters_frozen_before_entry;
  assert.deepEqual([p.a_short_lots, p.b_long_lots, p.a_multiplier, p.b_multiplier], [5, 8, 10, 5]);
  nearlyEqual(p.b_long_lots * p.b_multiplier, params.beta * p.a_short_lots * p.a_multiplier);
  assert.equal(s.entry.residual, s.entry.a - params.beta * s.entry.b);
  assert.equal(s.entry.z, 2);
  assert.deepEqual(s.scenarios.map(r => [r.short_a_pnl, r.long_b_pnl, r.net_pnl]), [[-2500, 4000, 1300], [-2000, 0, -2200]]);
  for (const r of s.scenarios) {
    assert.equal(r.short_a_pnl, (s.entry.a - r.a_exit) * p.a_short_lots * p.a_multiplier);
    assert.equal(r.long_b_pnl, (r.b_exit - s.entry.b) * p.b_long_lots * p.b_multiplier);
    assert.equal(r.short_a_pnl + r.long_b_pnl, r.gross_pnl);
    assert.equal(r.gross_pnl, p.short_residual_exposure * (s.entry.residual - r.exit_residual));
    assert.equal(r.gross_pnl - r.assumed_total_trading_cost, r.net_pnl);
    assert.equal(r.exit_z, (r.exit_residual - params.mean) / params.standard_deviation);
  }
  assert.equal(s.break_even_exit_residual, 96);
  assert.equal(p.short_residual_exposure * (s.entry.residual - s.break_even_exit_residual), s.scenarios[0].assumed_total_trading_cost);
});

test('a multiplier-matched statistical hedge cancels its assumed common component but not an arbitrary equal-price move', () => {
  const s = ledger.statistical_spread;
  const p = s.positions;
  const { beta } = s.parameters_frozen_before_entry;
  const matchedCommonMove = p.a_price_pnl_coefficient * beta * 100 + p.b_price_pnl_coefficient * 100;
  const equalPriceMove = p.a_price_pnl_coefficient * 100 + p.b_price_pnl_coefficient * 100;
  assert.equal(matchedCommonMove, 0);
  assert.equal(equalPriceMove, -1000);
  // Taking only 4 B lots from 0.8 × 5 ignores B's smaller multiplier.
  const wrongLotsCommonMove = p.a_price_pnl_coefficient * beta * 100 + 4 * p.b_multiplier * 100;
  assert.equal(wrongLotsCommonMove, -2000);
});

test('AR(1) half-life halves the conditional deviation rather than reaching the mean or guaranteeing an exit', () => {
  const s = ledger.statistical_spread;
  const { ar1_phi: phi, mean } = s.parameters_frozen_before_entry;
  nearlyEqual(phi ** s.half_life_in_bars, 0.5);
  nearlyEqual(s.conditional_mean_at_half_life - mean, (s.entry.residual - mean) / 2);
  assert.equal(s.conditional_mean_at_half_life, 80);
  assert.ok(s.half_life_in_bars > 6 && s.half_life_in_bars < 7);
  assert.ok(mean + phi ** 6 * (s.entry.residual - mean) > 80);
  assert.ok(mean + phi ** 7 * (s.entry.residual - mean) < 80);
  // The separate adverse scenario remains possible; the expectation is not a path constraint.
  assert.ok(s.scenarios[1].exit_residual > s.entry.residual);
});
