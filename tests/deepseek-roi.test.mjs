import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { defaults, evaluate, billTokens } from '../content/assets/deepseek-roi/model.js';
import { evaluateCashFlows } from '../content/assets/llm-economics/finance.js';
import { evaluateTokenScenario, billTokens as billExplicitTokens } from '../content/assets/llm-economics/token-model.js';
import { examplePreset, customDefaults } from '../content/assets/llm-economics/presets.js';

function near(actual, expected, relative = 1e-10) {
  assert.ok(Math.abs(actual - expected) <= relative * Math.max(1, Math.abs(expected)), `${actual} ≠ ${expected}`);
}
const highDemand = { monthlyInputTokens: 64e9, monthlyOutputTokens: 16e9 };
const demandAt = multiplier => ({ monthlyInputTokens: 24e9 * multiplier, monthlyOutputTokens: 6e9 * multiplier });
const annuity = (horizon = 36, rate = 0.1) => Array.from({ length: horizon }, (_, i) => Math.exp(-Math.log(1 + rate) * (i + 1) / 12)).reduce((sum, value) => sum + value, 0);
const recorded = JSON.parse(await readFile(new URL('../content/assets/deepseek-roi/scenarios.json', import.meta.url), 'utf8'));

test('API token categories and bucket tariffs include output exactly once', () => {
  const bill = billTokens(1e6, 1e6);
  assert.deepEqual(bill, { hit: { tokens: 750000, amount: 0.03 }, miss: { tokens: 250000, amount: 0.5 }, output: { tokens: 1e6, amount: 8 }, total: 8.53 });
  near(billTokens(1e6, 1e6, {}, 'off').total, 4.265);
  near(billTokens(1e6, 1e6, { apiCacheHit: 0 }).total, 10);
  near(billTokens(1e6, 1e6, { apiCacheHit: 1 }).total, 8.04);
  near(billTokens(0, 1e6).total, 8);
  near(billTokens(1e6, 0).total, 0.53);
  near(billTokens(1e6, 1e6, { prices: { peak: { output: 0 } } }).total, 0.53);
  for (const hit of [0, 1]) near(evaluate({ apiCacheHit: hit }).firstMonth.localInputTokens, 23.52e9);
});

test('fallback cache is independent and null follows the baseline rate', () => {
  for (const apiCacheHit of [0, 0.4, 1]) {
    near(billTokens(1e6, 1e6, { apiCacheHit }, 'peak', 1, 'fallback').total, billTokens(1e6, 1e6, { apiCacheHit }).total);
  }
  const independent = { apiCacheHit: 0.75, fallbackCacheHit: 0 };
  near(billTokens(1e6, 1e6, independent).total, 8.53);
  near(billTokens(1e6, 1e6, independent, 'peak', 1, 'fallback').total, 10);
  near(billTokens(1e6, 1e6, independent, 'off', 1, 'fallback').total, 5);
  near(billTokens(1e6, 1e6, { ...independent, apiAnnualChange: -0.3 }, 'peak', 13, 'fallback').total, 7);
});

test('cold fallback raises cost without changing baseline bills or capacity', () => {
  const warm = evaluate(highDemand);
  const cold = evaluate({ ...highDemand, fallbackCacheHit: 0 });
  near(cold.fullApiTco, warm.fullApiTco);
  near(cold.totalLocalInputTokens, warm.totalLocalInputTokens);
  const penalty = (18_091_212_800 * 1.96 + 640_000_000 * 0.98) * 0.75 / 1e6;
  near(cold.firstMonth.cachePenalty, penalty);
  near(cold.hybridTco - warm.hybridTco, penalty * 36);
  near(warm.NPV - cold.NPV, penalty * annuity());
  assert.ok(warm.NPV > 0 && cold.NPV < 0);
  for (const row of [...cold.months, ...cold.firstMonth.buckets]) {
    near(row.localDisplacedApiCost - row.cachePenalty, row.avoidedApiCost);
    near(row.fullApiCost - row.overflowApiCost, row.avoidedApiCost);
  }
});

test('zero-capacity fallback can have a cache penalty despite displacing no token bill', () => {
  const r = evaluate({ capacityMultiplier: 0, fallbackCacheHit: 0 });
  assert.equal(r.firstMonth.localDisplacedApiCost, 0);
  near(r.firstMonth.cachePenalty, (12e9 * 1.96 + 12e9 * 0.98) * 0.75 / 1e6);
  near(r.firstMonth.avoidedApiCost, -r.firstMonth.cachePenalty);
  assert.ok(evaluate({ fallbackCacheHit: 1 }).firstMonth.cachePenalty < 0);
});

test('zero workload has finite zero bills and token counts, but still incurs fixed cash costs', () => {
  const r = evaluate({ monthlyInputTokens: 0, monthlyOutputTokens: 0 });
  assert.equal(r.fullApiTco, 0);
  assert.equal(r.totalLocalInputTokens, 0);
  assert.equal(r.totalOverflowOutputTokens, 0);
  assert.equal(r.firstMonth.localInputShare, null);
  assert.equal(r.firstMonth.localOutputShare, null);
  assert.equal(r.capacity.maximumNoCongestionDemandMultiplier, null);
  assert.equal(r.simplePaybackMonth, null);
  near(r.hybridTco, 1_600_000 + 36 * 24_000 - 240_000);
  for (const b of r.firstMonth.buckets) assert.equal(b.localFraction, 0);
  assert.ok(!JSON.stringify(r).includes('NaN'));
});

test('zero-demand dimensions do not create zero-denominator or spurious capacity constraints', () => {
  const outputOnly = evaluate({ monthlyInputTokens: 0, inputTokensPerSecond: 0 });
  near(outputOnly.firstMonth.localOutputTokens, 5.88e9);
  assert.equal(outputOnly.firstMonth.localInputShare, null);
  near(outputOnly.firstMonth.baselineBillBreakdown.hit.tokens, 0);
  const inputOnly = evaluate({ monthlyOutputTokens: 0, outputTokensPerSecond: 0 });
  near(inputOnly.firstMonth.localInputTokens, 23.52e9);
  assert.equal(inputOnly.firstMonth.localOutputShare, null);
  assert.equal(inputOnly.firstMonth.baselineBillBreakdown.output.amount, 0);
});

test('zero capacity in a demanded dimension sends both coupled token streams to API', () => {
  for (const override of [{ capacityMultiplier: 0 }, { inputTokensPerSecond: 0 }, { outputTokensPerSecond: 0 }, { availability: 0 }, { headroom: 0 }]) {
    const r = evaluate(override);
    assert.equal(r.totalLocalInputTokens, 0);
    assert.equal(r.totalLocalOutputTokens, 0);
    near(r.totalOverflowInputTokens, 36 * 24e9);
    near(r.totalOverflowOutputTokens, 36 * 6e9);
    near(r.totalOverflowApiCost, r.fullApiTco);
    assert.equal(r.simplePaybackMonth, null);
  }
});

test('input and output can constrain different buckets while each bucket routes proportionally', () => {
  const r = evaluate({ monthlyInputTokens: 1e6, monthlyOutputTokens: 1e6, inputTokensPerSecond: 100, outputTokensPerSecond: 200, peakHours: 1, offHours: 1, headroom: 1, availability: 0.8, inputPeakShare: 0.8, outputPeakShare: 0.2 });
  const [peak, off] = r.firstMonth.buckets;
  near(peak.onlineFraction, 0.45); // Input: 360k / 800k.
  near(off.onlineFraction, 0.9); // Output: 720k / 800k.
  near(peak.localInputTokens, 288000);
  near(peak.localOutputTokens, 72000);
  near(off.localInputTokens, 144000);
  near(off.localOutputTokens, 576000);
  near(r.firstMonth.localInputShare, 0.432);
  near(r.firstMonth.localOutputShare, 0.648);
  for (const b of [peak, off]) near(b.localInputTokens / b.demandInputTokens, b.localOutputTokens / b.demandOutputTokens);
  // More input-only throughput cannot help the output-constrained off bucket.
  near(evaluate({ ...r.parameters, inputTokensPerSecond: 200 }).firstMonth.buckets[1].localOutputTokens, off.localOutputTokens);
});

test('buckets do not lend off-peak spare capacity to peak traffic and conserve both token streams', () => {
  const r = evaluate({ monthlyInputTokens: 48e9, monthlyOutputTokens: 12e9 });
  const [peak, off] = r.firstMonth.buckets;
  near(peak.effectiveInputTokenCapacity, 154 * 3600 * 32000 * 0.8 * 0.98);
  near(peak.localInputTokens, 13_908_787_200);
  near(peak.localOutputTokens, 3_477_196_800);
  near(off.overflowInputTokens, 480_000_000);
  near(off.downtimeOverflowInputTokens, 480_000_000);
  assert.equal(off.overloadOverflowInputTokens, 0);
  assert.ok(off.localInputTokens < off.effectiveInputTokenCapacity);
  for (const row of [...r.months, ...r.firstMonth.buckets]) {
    for (const stream of ['Input', 'Output']) {
      near(row[`local${stream}Tokens`] + row[`overflow${stream}Tokens`], row[`demand${stream}Tokens`]);
      near(row[`downtimeOverflow${stream}Tokens`] + row[`overloadOverflow${stream}Tokens`], row[`overflow${stream}Tokens`]);
    }
    for (const bill of [row.baselineBillBreakdown, row.fallbackBillBreakdown]) near(bill.total, bill.hit.amount + bill.miss.amount + bill.output.amount);
    near(row.baselineBillBreakdown.hit.tokens + row.baselineBillBreakdown.miss.tokens, row.demandInputTokens);
    near(row.fallbackBillBreakdown.output.tokens, row.overflowOutputTokens);
  }
  near(r.totalDemandInputTokens, r.totalLocalInputTokens + r.totalOverflowInputTokens);
  near(r.totalDemandOutputTokens, r.totalLocalOutputTokens + r.totalOverflowOutputTokens);
});

test('availability is applied after online congestion; low traffic still has downtime fallback', () => {
  const r = evaluate({ monthlyInputTokens: 4e9, monthlyOutputTokens: 1e9 });
  near(r.firstMonth.localInputTokens, 3.92e9);
  near(r.firstMonth.overflowInputTokens, 80e6);
  near(r.firstMonth.overflowOutputTokens, 20e6);
  near(r.firstMonth.overloadOverflowInputTokens, 0);
  assert.equal(r.capacity.maximumAllLocalDemandMultiplier, 0);
  const alwaysOnline = evaluate({ ...r.parameters, availability: 1 });
  assert.equal(alwaysOnline.firstMonth.overflowInputTokens, 0);
  near(alwaysOnline.capacity.maximumNoCongestionDemandMultiplier, r.capacity.maximumNoCongestionDemandMultiplier);
  near(alwaysOnline.capacity.maximumAllLocalDemandMultiplier, r.capacity.maximumNoCongestionDemandMultiplier);
  const multiplier = evaluate().capacity.maximumNoCongestionDemandMultiplier;
  near(multiplier, 1.18272);
  const atBoundary = evaluate(demandAt(multiplier));
  near(atBoundary.firstMonth.overloadOverflowInputTokens, 0);
  near(atBoundary.firstMonth.downtimeOverflowInputTokens, 24e9 * multiplier * 0.02);
  const above = evaluate(demandAt(multiplier + 1e-6));
  assert.ok(above.firstMonth.overloadOverflowInputTokens > 0);
  near(above.firstMonth.buckets[0].localInputTokens, 13_908_787_200);
});

test('with fixed mix and equal cache rates, throughput saturation caps savings', () => {
  const saturated = evaluate(demandAt(5));
  const doubled = evaluate(demandAt(10));
  near(saturated.firstMonth.localInputTokens, 65_028_096_000);
  near(saturated.firstMonth.localOutputTokens, 16_257_024_000);
  near(saturated.firstMonth.operatingSavings, 75_855.157248);
  near(saturated.NPV, 945_969.26135909);
  near(saturated.NPV, doubled.NPV);
  const cold = evaluate({ ...demandAt(5), fallbackCacheHit: 0 });
  assert.ok(evaluate({ ...demandAt(10), fallbackCacheHit: 0 }).NPV < cold.NPV);
});

test('default and fixed token-volume cases preserve independently calculated financial results', () => {
  const r = evaluate();
  near(r.NPV, -1_600_000 + (45_540 * 0.98 - 24_000) * annuity() + 240_000 / 1.1 ** 3);
  near(r.fullApiTco, 45_540 * 36);
  near(r.hybridTco, 1_600_000 + (24_000 + 45_540 * 0.02) * 36 - 240_000);
  near(r.terminalNetCash, -617_348.8);
  const golden = [
    [{}, -776332.759287404, 2256788.8],
    [{ monthlyInputTokens: 4e9, monthlyOutputTokens: 1e9 }, -1936188.3352119583, 2229464.8],
    [{ monthlyInputTokens: 48e9, monthlyOutputTokens: 12e9 }, -142847.48025123496, 3164966.861824],
    [highDemand, 166447.339995313, 3900893.261824],
    [{ ...highDemand, fallbackCacheHit: 0 }, -677598.0286060113, 4875214.6432],
    [{ apiAnnualChange: -0.3 }, -1279042.3204302639, 2244433.1488042474],
  ];
  for (const [overrides, npv, hybridTco] of golden) {
    const result = evaluate(overrides);
    near(result.NPV, npv);
    near(result.hybridTco, hybridTco);
  }
});

test('cash-flow identities hold with independent demand growth and extra token costs', () => {
  for (const overrides of [{}, highDemand, { apiAnnualChange: -0.3 }, { inputAnnualGrowth: 0.2, outputAnnualGrowth: -0.1 }, { fallbackCacheHit: 0, perMillionLocalInputCost: 0.1, perMillionLocalOutputCost: 0.4, opexAnnualGrowth: 0.2 }]) {
    const r = evaluate(overrides);
    near(r.fullApiTco - r.hybridTco, r.terminalNetCash);
    near(r.capitalROI * r.parameters.capex, r.terminalNetCash);
    near(r.months.reduce((sum, m) => sum + m.cashFlow, 0), r.terminalNetCash);
    near(r.months.reduce((sum, m) => sum + m.discountedCashFlow, 0), r.NPV);
    near(r.totalOpex, r.totalFixedOpex + r.totalVariableLocalCost);
    near(r.firstMonth.capitalAllocatedCost, (r.parameters.capex - r.parameters.terminalResidual) / r.parameters.horizon);
    near(r.firstMonth.grossAllocatedHybridCost, r.firstMonth.hybridOperatingCost + r.parameters.capex / r.parameters.horizon);
  }
  assert.equal(evaluate({ capex: 0 }).capitalROI, null);
});

test('extra per-million token costs apply once only to locally served tokens', () => {
  const base = evaluate();
  const r = evaluate({ perMillionLocalInputCost: 0.1, perMillionLocalOutputCost: 0.1 });
  near(r.firstMonth.variableLocalInputCost, 2352);
  near(r.firstMonth.variableLocalOutputCost, 588);
  near(r.firstMonth.variableLocalCost, 2940);
  near(r.firstMonth.fixedOpex, 24000);
  near(r.firstMonth.opex, 26940);
  near(r.hybridTco - base.hybridTco, 2940 * 36);
  near(r.totalOverflowApiCost, base.totalOverflowApiCost);
  assert.equal(evaluate({ capacityMultiplier: 0, perMillionLocalOutputCost: 1 }).totalVariableLocalCost, 0);
  assert.equal(evaluate({ monthlyInputTokens: 0, monthlyOutputTokens: 0, perMillionLocalInputCost: 1 }).totalVariableLocalCost, 0);
});

test('prices, input/output growth and fixed/token OpEx growth compound independently by month', () => {
  const r = evaluate({ monthlyInputTokens: 1e6, monthlyOutputTokens: 2e6, inputPeakShare: 0.25, outputPeakShare: 0.75, inputAnnualGrowth: 0.2, outputAnnualGrowth: -0.2, apiAnnualChange: -0.3, monthlyOpex: 100, opexAnnualGrowth: 0.1, perMillionLocalInputCost: 1, perMillionLocalOutputCost: 2 });
  near(r.months[1].fullApiCost, 0.33125 + 14);
  near(r.months[13].demandInputTokens, 1.2e6);
  near(r.months[13].demandOutputTokens, 1.6e6);
  near(r.months[13].fullApiCost, (0.33125 * 1.2 + 14 * 0.8) * 0.7);
  near(r.months[13].fixedOpex, 110);
  near(r.months[13].variableLocalCost, 0.98 * (1.2 + 1.6 * 2) * 1.1);
  near(r.months[2].demandInputTokens, 1e6 * 1.2 ** (1 / 12));
  near(r.monthlyDiscount, 1.1 ** (1 / 12) - 1);
  near(r.months[1].discountFactor, 1.1 ** (-1 / 12));
  near(r.months[12].discountFactor, 1 / 1.1);
});

test('arbitrary horizons and custom tariffs determine cash flows without 36-month constants', () => {
  for (const horizon of [1, 12, 60]) {
    const r = evaluate({ horizon, capex: 100, monthlyOpex: 2, terminalResidual: 10, monthlyInputTokens: 1e6, monthlyOutputTokens: 2e6, availability: 1, prices: { peak: { hit: 3, miss: 5, output: 7 }, off: { hit: 2, miss: 4, output: 6 } } });
    const fullApiCost = 0.5 * (0.75 * 3 + 0.25 * 5) + 0.5 * (0.75 * 2 + 0.25 * 4) + 7 + 6;
    near(fullApiCost, 16);
    near(r.NPV, -100 + (16 - 2) * annuity(horizon) + 10 * 1.1 ** (-horizon / 12));
    near(r.fullApiTco, 16 * horizon);
    assert.equal(r.months.length, horizon + 1);
    near(r.months.at(-1).residual, 10);
  }
  const zeroRate = evaluate({ annualDiscount: 0 });
  near(zeroRate.NPV, zeroRate.terminalNetCash);
});

test('terminal proceeds appear only at the chosen horizon and never create operating payback', () => {
  const r = evaluate({ terminalResidual: 2e6, horizon: 24 });
  assert.ok(r.terminalNetCash > 0);
  assert.equal(r.simplePaybackMonth, null);
  assert.equal(r.discountedPaybackMonth, null);
  for (const row of r.months.slice(0, -1)) {
    assert.equal(row.residual, 0);
    near(row.cumulativeNet, row.cumulativeOperatingNet);
  }
  near(r.terminalNetCash - r.operatingNetCash, 2e6);
  near(r.NPV - r.operatingNPV, 2e6 / 1.1 ** 2);
});

test('negative savings never repay, and an early crossing followed by losses is exposed', () => {
  const low = evaluate({ monthlyInputTokens: 4e9, monthlyOutputTokens: 1e9 });
  assert.ok(low.months.slice(1).every(m => m.operatingSavings < 0));
  assert.equal(low.simplePaybackMonth, null);
  assert.equal(low.discountedPaybackMonth, null);
  const reversed = evaluate({ capex: 10000, terminalResidual: 0, apiAnnualChange: -0.9 });
  assert.equal(reversed.simplePaybackMonth, 1);
  assert.equal(reversed.discountedPaybackMonth, 1);
  assert.ok(reversed.simplePaybackReversalMonths.length > 0);
  assert.ok(reversed.discountedPaybackReversalMonths.length > 0);
  assert.equal(reversed.sustainedSimplePaybackMonth, null);
  assert.equal(reversed.sustainedDiscountedPaybackMonth, null);
});

test('payback is an integer month-end crossing and excludes final resale', () => {
  const high = evaluate(highDemand);
  assert.equal(high.simplePaybackMonth, 32);
  assert.equal(high.discountedPaybackMonth, null);
  const saturated = evaluate(demandAt(5));
  assert.equal(saturated.simplePaybackMonth, 22);
  assert.equal(saturated.discountedPaybackMonth, 24);
  assert.ok(saturated.months[21].cumulativeOperatingNet < 0);
  assert.ok(saturated.months[22].cumulativeOperatingNet > 0);
});

test('investment limits and hypothetical service ceiling retain their financial meanings', () => {
  for (const overrides of [{}, { horizon: 12 }, { ...highDemand, fallbackCacheHit: 0 }, { apiAnnualChange: -0.3 }]) {
    const r = evaluate(overrides);
    near(r.maxInitialInvestment - r.parameters.capex, r.NPV);
    near(r.maxInitialInvestmentWithoutResale, r.presentValueOperatingSavings);
    near(r.maxInitialInvestment - r.maxInitialInvestmentWithoutResale, r.presentValueResidual);
    near(r.idealServiceCeiling, r.presentValueFullApi + r.presentValueResidual);
    assert.ok(r.maxInitialInvestment <= r.idealServiceCeiling);
    if (r.maxInitialInvestment >= 0) near(evaluate({ ...overrides, capex: r.maxInitialInvestment }).NPV, 0, 1e-8);
  }
  const r = evaluate();
  near(r.maxInitialInvestment, 823667.240712596);
  near(evaluate({ terminalResidual: 0 }).maxInitialInvestment, r.maxInitialInvestmentWithoutResale);
  const proportional = r.maxInitialInvestmentWithoutResale / (1 - 0.15 / 1.1 ** 3);
  near(proportional, 725064.4347065752);
  near(evaluate({ capex: proportional, terminalResidual: 0.15 * proportional }).NPV, 0, 1e-8);
  near(evaluate({ availability: 1, monthlyOpex: 0, capacityMultiplier: 10 }).maxInitialInvestment, r.idealServiceCeiling);
  near(evaluate({ monthlyOpex: 0 }).NPV, -27857.75667415664);
});

test('capacity expansion gains depend on an actual binding token constraint', () => {
  near(evaluate({ capacityMultiplier: 2 }).NPV, evaluate().NPV);
  const mid = { monthlyInputTokens: 48e9, monthlyOutputTokens: 12e9 };
  near(evaluate({ ...mid, capacityMultiplier: 1.25 }).firstMonth.operatingSavings - evaluate(mid).firstMonth.operatingSavings, 8797.307904);
  near(evaluate({ ...mid, capacityMultiplier: 1.25 }).NPV - evaluate(mid).NPV, 274356.877351498);
});

test('schema5 assets exactly reproduce evaluations and token-volume thresholds straddle zero', () => {
  assert.equal(recorded.schemaVersion, 5);
  for (const scenario of recorded.scenarios) assert.deepEqual(scenario.result, evaluate(scenario.overrides));
  for (const [name, metric] of [['npvBreakEvenIncludingTerminal', 'NPV'], ['operatingPaybackWithinHorizon', 'operatingNetCash'], ['discountedOperatingPaybackWithinHorizon', 'operatingNPV']]) {
    const boundary = recorded.thresholds[name];
    near(boundary.monthlyInputTokens, boundary.demandMultiplier * 24e9);
    near(boundary.monthlyOutputTokens, boundary.demandMultiplier * 6e9);
    assert.ok(evaluate(demandAt(boundary.demandMultiplier - 1e-6))[metric] < 0);
    assert.ok(evaluate(demandAt(boundary.demandMultiplier + 1e-6))[metric] > 0);
  }
  near(recorded.independentChecks.analyticNpvThreshold, recorded.thresholds.npvBreakEvenIncludingTerminal.demandMultiplier);
  assert.equal(recorded.independentChecks.passed, true);
});

test('rental memo replaces ownership OpEx, separately prices fallback, and breaks even at Rmax', () => {
  const memo = recorded.decisionMemo.rental;
  for (const c of [memo.base, memo.high]) {
    const demand = { monthlyInputTokens: c.monthlyInputTokens, monthlyOutputTokens: c.monthlyOutputTokens };
    const inputs = { ...demand, capex: 100000, terminalResidual: 0 };
    near(evaluate({ ...inputs, monthlyOpex: c.maxMonthlyRentalRate }).NPV, 0, 1e-8);
    near(evaluate({ ...inputs, monthlyOpex: c.maxMonthlyRentalRate + 1 }).NPV, -annuity(), 1e-7);
    const testRate = 60000; // Arithmetic probe, not a supplier quote.
    const lease = evaluate({ ...inputs, monthlyOpex: testRate });
    near(lease.hybridTco, c.rentalTcoIntercept + 36 * testRate);
    near(lease.NPV, c.rentalNPVIntercept - annuity() * testRate);
    near(lease.totalOpex, 36 * testRate);
    near(c.maxMonthlyRentalRateZeroSetup, c.monthlyAvoidedApiCost);
    near(evaluate({ ...inputs, monthlyOpex: memo.monthlyRateEqualToOwnershipPresentValue }).NPV, evaluate(demand).NPV);
  }
  near(memo.base.maxMonthlyRentalRate, 41422.680087350185);
  near(memo.high.maxMonthlyRentalRate, 71653.11170335018);
});

test('high-demand fragility boundaries match independent token arithmetic and have opposite-sign probes', () => {
  const d = recorded.decisionMemo.fragilityHigh;
  const required = 24000 + (1600000 - 240000 / 1.1 ** 3) / annuity();
  const peakAvoided = (13_908_787_200 * 0.53 + 3_477_196_800 * 8) / 1e6;
  const offAvoided = (31.36e9 * 0.265 + 7.84e9 * 4) / 1e6;
  near(d.immediateApiPriceReduction.value, 1 - required / (peakAvoided + offAvoided));
  near(d.minimumCapacityMultiplier.value, (required - offAvoided) / peakAvoided);
  near(d.minimumCapacityMultiplier.inputTokensPerSecond, d.minimumCapacityMultiplier.value * 32000);
  near(d.minimumCapacityMultiplier.outputTokensPerSecond, d.minimumCapacityMultiplier.value * 8000);
  const penaltyPerFraction = (18_091_212_800 * 1.96 + 640_000_000 * 0.98) / 1e6;
  const highNPV = (peakAvoided + offAvoided - 24000) * annuity() - 1600000 + 240000 / 1.1 ** 3;
  near(d.minimumFallbackCacheHit.value, 0.75 - highNPV / (penaltyPerFraction * annuity()));
  const pricesAt = reduction => Object.fromEntries(Object.entries(defaults.prices).map(([bucket, prices]) => [bucket, Object.fromEntries(Object.entries(prices).map(([type, price]) => [type, price * (1 - reduction)]))]));
  for (const [boundary, fn, increasing] of [
    [d.immediateApiPriceReduction, x => evaluate({ ...highDemand, prices: pricesAt(x) }).NPV, false],
    [d.minimumCapacityMultiplier, x => evaluate({ ...highDemand, capacityMultiplier: x }).NPV, true],
    [d.minimumFallbackCacheHit, x => evaluate({ ...highDemand, fallbackCacheHit: x }).NPV, true],
  ]) {
    near(fn(boundary.value), 0, 1e-8);
    near(fn(boundary.belowValue), boundary.NPVbelow, 1e-8);
    near(fn(boundary.aboveValue), boundary.NPVabove, 1e-8);
    assert.ok(increasing ? boundary.NPVbelow < 0 && boundary.NPVabove > 0 : boundary.NPVbelow > 0 && boundary.NPVabove < 0);
  }
});

test('high-demand stress cases distinguish annual price decline from a one-time token-demand reduction', () => {
  const cases = recorded.decisionMemo.stressesHigh;
  for (const c of [cases.apiAnnualMinus30, cases.monthlyDemandMinus20, cases.combined]) {
    const r = evaluate(c.overrides);
    near(r.NPV, c.NPV);
    near(r.fullApiTco, c.fullApiTco);
    near(r.hybridTco, c.hybridTco);
    assert.ok(r.NPV < 0);
  }
  const reduced = evaluate(cases.monthlyDemandMinus20.overrides);
  assert.ok(reduced.months.slice(1).every(m => m.demandInputTokens === 51.2e9 && m.demandOutputTokens === 12.8e9));
  const annual = evaluate(cases.apiAnnualMinus30.overrides);
  near(annual.months[13].avoidedApiCost / annual.months[1].avoidedApiCost, 0.7);
  near(cases.apiAnnualMinus30.NPV, -676782.0380127078);
  near(cases.monthlyDemandMinus20.NPV, -80988.51620192558);
  near(cases.combined.NPV, -834847.3055623267);
});

test('invalid and obsolete interfaces fail explicitly rather than silently preserving a hidden volume model', () => {
  for (const overrides of [{ monthlyRequests: 1 }, { inputTokens: 1 }, { outputTokens: 1 }, { qSlo: 1 }, { perLocalRequestCost: 1 }, { monthlyInputTokens: -1 }, { outputTokensPerSecond: NaN }, { capacityMultiplier: -1 }, { apiCacheHit: 2 }, { fallbackCacheHit: undefined }, { perMillionLocalInputCost: -1 }, { perMillionLocalOutputCost: -1 }, { annualDiscount: -1 }, { inputAnnualGrowth: -1 }, { outputAnnualGrowth: -1 }, { opexAnnualGrowth: -1 }, { horizon: 1.5 }, { peakHours: 0, offHours: 0 }, { prices: { peak: { output: -1 } } }]) assert.throws(() => evaluate(overrides), RangeError);
  assert.throws(() => billTokens(-1, 0), RangeError);
  assert.throws(() => billTokens(0, 1, {}, 'night'), RangeError);
  assert.throws(() => billTokens(0, 1, {}, 'peak', 0), RangeError);
  assert.throws(() => billTokens(0, 1, {}, 'peak', 1, 'local'), RangeError);
});

test('generic finance matches a hand-discounted series with later capital expenditure and final proceeds', () => {
  // Synthetic currency units; month 2 includes a replacement outlay.
  // The annual rate is chosen to produce exactly 10% per month for easy hand checks.
  const cash = evaluateCashFlows({
    initialInvestment: 100, monthlyNetCashFlows: [60, -40, 100],
    terminalValue: 20, annualDiscount: 1.1 ** 12 - 1,
  });
  near(cash.monthlyDiscount, 0.1);
  near(cash.preTerminalNetCash, 20);
  near(cash.terminalNetCash, 40);
  near(cash.capitalROI, 0.4);
  near(cash.NPV, -100 + 60 / 1.1 - 40 / 1.1 ** 2 + 120 / 1.1 ** 3);
  near(cash.preTerminalNPV, -100 + 60 / 1.1 - 40 / 1.1 ** 2 + 100 / 1.1 ** 3);
  near(cash.presentValueTerminalValue, 20 / 1.1 ** 3);
  near(cash.maxInitialInvestment - 100, cash.NPV);
  assert.equal(cash.simplePaybackMonth, 3);
  assert.equal(cash.discountedPaybackMonth, null);
  assert.equal(cash.months[2].terminalValue, 0);
  assert.equal(cash.months[3].terminalValue, 20);
});

test('generic finance exposes payback reversal and permits explicit disposal costs', () => {
  const inputs = { initialInvestment: 100, monthlyNetCashFlows: [120, -80, 100], terminalValue: -10 };
  const cash = evaluateCashFlows(inputs);
  near(cash.NPV, 30);
  near(cash.preTerminalNetCash, 40);
  assert.equal(cash.simplePaybackMonth, 1);
  assert.deepEqual(cash.simplePaybackReversalMonths, [2]);
  assert.equal(cash.sustainedSimplePaybackMonth, 3);
  assert.deepEqual(cash.discountedPaybackReversalMonths, [2]);
  near(cash.months.reduce((sum, row) => sum + row.cashFlow, 0), cash.terminalNetCash);
  const noInvestment = evaluateCashFlows({ initialInvestment: 0, monthlyNetCashFlows: [1] });
  assert.equal(noInvestment.capitalROI, null);
  assert.equal(noInvestment.simplePaybackMonth, 0);
  assert.deepEqual(inputs.monthlyNetCashFlows, [120, -80, 100]);
});

test('generic cash-flow validation rejects missing, nonfinite or invalid financial inputs', () => {
  for (const parameters of [
    undefined, {}, { initialInvestment: 1, monthlyNetCashFlows: [] },
    { initialInvestment: -1, monthlyNetCashFlows: [1] },
    { initialInvestment: 1, monthlyNetCashFlows: [NaN] },
    { initialInvestment: 1, monthlyNetCashFlows: [Infinity] },
    { initialInvestment: 1, monthlyNetCashFlows: [1], terminalValue: Infinity },
    { initialInvestment: 1, monthlyNetCashFlows: [1], annualDiscount: -1 },
  ]) assert.throws(() => evaluateCashFlows(parameters), RangeError);
});

// All values here are invented arithmetic inputs, unrelated to a model or vendor.
const syntheticParameters = {
  monthlyInputTokens: 1e6, monthlyOutputTokens: 0.5e6,
  inputTokensPerSecond: 200, outputTokensPerSecond: 100, capacityMultiplier: 1,
  inputPeakShare: 0.75, outputPeakShare: 0.5,
  inputAnnualGrowth: 0, outputAnnualGrowth: 0, opexAnnualGrowth: 0,
  perMillionLocalInputCost: 0.2, perMillionLocalOutputCost: 0.4,
  apiCacheHit: 0.5, fallbackCacheHit: 0.25,
  peakHours: 1, offHours: 1, availability: 0.8, headroom: 0.5,
  capex: 10, monthlyOpex: 0.1, horizon: 2, annualDiscount: 0, terminalResidual: 1,
  apiAnnualChange: 0,
  prices: { peak: { hit: 1, miss: 3, output: 4 }, off: { hit: 0.5, miss: 1.5, output: 2 } },
};

test('an explicit non-vendor token scenario matches hand billing, routing and cash arithmetic', () => {
  const r = evaluateTokenScenario(syntheticParameters);
  const [peak, off] = r.firstMonth.buckets;
  near(peak.localFraction, 0.384); // Input constraint: .8 × (360k / 750k).
  near(off.localFraction, 0.576); // Output constraint: .8 × (180k / 250k).
  near(r.firstMonth.localInputTokens, 432000);
  near(r.firstMonth.localOutputTokens, 240000);
  near(r.firstMonth.fullApiCost, 2.5 + 0.75);
  near(r.firstMonth.overflowApiCost, 1.771 + 0.3445);
  near(r.firstMonth.localDisplacedApiCost, 1.392);
  near(r.firstMonth.cachePenalty, 0.2575);
  near(r.firstMonth.variableLocalCost, 0.432 * 0.2 + 0.240 * 0.4);
  near(r.firstMonth.operatingSavings, 0.8521);
  near(r.fullApiTco, 6.5);
  near(r.hybridTco, 13.7958);
  near(r.NPV, -7.2958);
  near(r.capitalROI, -0.72958);
  near(r.maxInitialInvestment, 2.7042);
  assert.equal(r.simplePaybackMonth, null);
  near(billExplicitTokens(0.75e6, 0.25e6, syntheticParameters).total, 2.5);
  const finance = evaluateCashFlows({ initialInvestment: 10, monthlyNetCashFlows: [0.8521, 0.8521], terminalValue: 1 });
  near(r.NPV, finance.NPV);
  near(r.capitalROI, finance.capitalROI);
});

test('generic token calculation requires a full explicit profile and never chooses the example', () => {
  assert.throws(() => evaluateTokenScenario(), /explicit token profile/);
  assert.throws(() => evaluateTokenScenario({}), /Missing token-model parameter/);
  const { monthlyOpex, ...missingCost } = syntheticParameters;
  assert.throws(() => evaluateTokenScenario(missingCost), /Missing token-model parameter: monthlyOpex/);
  assert.throws(() => evaluateTokenScenario({ ...syntheticParameters, prices: { peak: { hit: 1 } } }), /Missing prices/);
  assert.throws(() => evaluateTokenScenario({ ...syntheticParameters, metadata: { currency: 'USD' } }), /Unknown token-model parameter/);
  const blank = evaluateTokenScenario(customDefaults);
  near(blank.NPV, 0);
  near(blank.fullApiTco, 0);
  near(blank.hybridTco, 0);
  assert.equal(blank.capitalROI, null);
  for (const b of ['peak', 'off']) {
    assert.deepEqual(customDefaults.prices[b], { hit: 0, miss: 0, output: 0 });
  }
  assert.deepEqual(evaluateTokenScenario(examplePreset.parameters), evaluate());
});

test('model labels have no numerical effect and all monetary inputs share one scalable unit', () => {
  const one = { metadata: { name: 'Synthetic A', version: '1', currency: 'USD' }, parameters: syntheticParameters };
  const two = { metadata: { name: 'Synthetic B', version: '2', currency: 'CNY' }, parameters: syntheticParameters };
  // Relabeling is not currency conversion and cannot change the numeric result.
  assert.deepEqual(evaluateTokenScenario(one.parameters), evaluateTokenScenario(two.parameters));
  const factor = 7;
  const scaled = { ...syntheticParameters };
  for (const key of ['capex', 'monthlyOpex', 'terminalResidual', 'perMillionLocalInputCost', 'perMillionLocalOutputCost']) scaled[key] *= factor;
  scaled.prices = Object.fromEntries(Object.entries(syntheticParameters.prices).map(([bucket, prices]) => [bucket, Object.fromEntries(Object.entries(prices).map(([type, price]) => [type, price * factor]))]));
  const a = evaluateTokenScenario(syntheticParameters);
  const b = evaluateTokenScenario(scaled);
  for (const key of ['NPV', 'fullApiTco', 'hybridTco', 'maxInitialInvestment', 'presentValueFullApi']) near(b[key], a[key] * factor);
  near(b.capitalROI, a.capitalROI);
  assert.equal(b.simplePaybackMonth, a.simplePaybackMonth);
  near(b.firstMonth.localInputTokens, a.firstMonth.localInputTokens);
});
