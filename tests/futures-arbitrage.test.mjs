import test from 'node:test';
import assert from 'node:assert/strict';
import { calculateCalendarSpread } from '../content/assets/futures-arbitrage/calculator.js';

const example = { nearOpen: 4000, farOpen: 4100, nearClose: 4050, farClose: 4110 };

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
