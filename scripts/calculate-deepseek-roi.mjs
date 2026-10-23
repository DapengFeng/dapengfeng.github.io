#!/usr/bin/env node
import assert from 'node:assert/strict';
import { readFile, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { defaults, evaluate, billTokens } from '../content/assets/deepseek-roi/model.js';
import { examplePreset } from '../content/assets/llm-economics/presets.js';
import { evaluateTokenScenario } from '../content/assets/llm-economics/token-model.js';

const destination = fileURLToPath(new URL('../content/assets/deepseek-roi/scenarios.json', import.meta.url));
const demandAt = demandMultiplier => ({ monthlyInputTokens: defaults.monthlyInputTokens * demandMultiplier, monthlyOutputTokens: defaults.monthlyOutputTokens * demandMultiplier });
const thresholdAt = demandMultiplier => demandMultiplier === null ? null : ({ demandMultiplier, ...demandAt(demandMultiplier) });
const midDemand = { monthlyInputTokens: 48e9, monthlyOutputTokens: 12e9 };
const highDemand = { monthlyInputTokens: 64e9, monthlyOutputTokens: 16e9 };
const definitions = [
  { id: 'base-24b-6b', label: { en: '24B input / 6B output tokens per month', zh: '每月 240 亿输入 / 60 亿输出 token' }, overrides: {} },
  { id: 'low-4b-1b', label: { en: '4B input / 1B output tokens per month', zh: '每月 40 亿输入 / 10 亿输出 token' }, overrides: { monthlyInputTokens: 4e9, monthlyOutputTokens: 1e9 } },
  { id: 'mid-48b-12b', label: { en: '48B input / 12B output tokens per month', zh: '每月 480 亿输入 / 120 亿输出 token' }, overrides: midDemand },
  { id: 'high-64b-16b', label: { en: '64B input / 16B output tokens per month', zh: '每月 640 亿输入 / 160 亿输出 token' }, overrides: highDemand },
  { id: 'high-fallback-cold', label: { en: '64B / 16B tokens; fallback cache hit 0%', zh: '每月 640 亿 / 160 亿 token，回退缓存命中率 0%' }, overrides: { ...highDemand, fallbackCacheHit: 0 } },
  { id: 'base-api-minus-30', label: { en: 'Base demand; API price −30%/year', zh: '基准 token 需求；API 每年降价 30%' }, overrides: { apiAnnualChange: -0.3 } },
  { id: 'half-capacity', label: { en: 'Simultaneous input/output capacity ×0.5', zh: '输入、输出同时可达容量均降至 0.5 倍' }, optional: true, overrides: { capacityMultiplier: 0.5 } },
  { id: 'capex-minus-30', label: { en: 'Initial investment −30%', zh: '初始投资降低 30%' }, optional: true, overrides: { capex: defaults.capex * 0.7 } },
];
function near(actual, expected, label) {
  assert.ok(Math.abs(actual - expected) <= 1e-8 * Math.max(1, Math.abs(expected)), `${label}: ${actual} ≠ ${expected}`);
}
function bisectIncreasing(fn, low, high) {
  assert.ok(fn(low) < 0 && fn(high) > 0, 'Boundary must be bracketed');
  for (let i = 0; i < 70; i++) {
    const midpoint = (low + high) / 2;
    if (fn(midpoint) >= 0) high = midpoint;
    else low = midpoint;
  }
  return (low + high) / 2;
}
const baseline = evaluate();
assert.deepEqual(baseline, evaluateTokenScenario(examplePreset.parameters), 'The compatibility entry point must use the explicit generic token profile');
const annuityFactor = Array.from({ length: defaults.horizon }, (_, i) => (1 + defaults.annualDiscount) ** (-(i + 1) / 12)).reduce((a, b) => a + b, 0);
const terminalDiscountFactor = (1 + defaults.annualDiscount) ** (-defaults.horizon / 12);
const residualPV = defaults.terminalResidual * terminalDiscountFactor;
const saturationDemandMultiplier = Math.max(...baseline.capacity.buckets.map(b => b.noCongestionDemandMultiplier ?? 0));
function demandThreshold(metric) {
  if (metric(evaluate(demandAt(0))) >= 0) return thresholdAt(0);
  if (metric(evaluate(demandAt(saturationDemandMultiplier))) < 0) return null;
  return thresholdAt(bisectIncreasing(m => metric(evaluate(demandAt(m))), 0, saturationDemandMultiplier));
}
const uncongestedAvoidedApiCost = 45_540 * defaults.availability;
const thresholds = {
  referenceDemand: demandAt(1),
  unconstrainedMonthlyOperatingBreakEven: thresholdAt(defaults.monthlyOpex / uncongestedAvoidedApiCost),
  unconstrainedOperatingPaybackWithinHorizon: thresholdAt((defaults.monthlyOpex + defaults.capex / defaults.horizon) / uncongestedAvoidedApiCost),
  monthlyOperatingBreakEven: demandThreshold(r => r.firstMonth.operatingSavings),
  operatingPaybackWithinHorizon: demandThreshold(r => r.operatingNetCash),
  discountedOperatingPaybackWithinHorizon: demandThreshold(r => r.operatingNPV),
  npvBreakEvenIncludingTerminal: demandThreshold(r => r.NPV),
  maximumAllLocal: thresholdAt(baseline.capacity.maximumAllLocalDemandMultiplier),
  maximumNoCongestion: thresholdAt(baseline.capacity.maximumNoCongestionDemandMultiplier),
  totalMonthlyInputTokenCapacity: baseline.capacity.totalMonthlyInputTokenCapacity,
  totalMonthlyOutputTokenCapacity: baseline.capacity.totalMonthlyOutputTokenCapacity,
  saturationByBucket: baseline.capacity.buckets.map(b => ({ id: b.id, ...thresholdAt(b.noCongestionDemandMultiplier) })),
};
// Independent token arithmetic: input-cache hit/miss and output are distinct bill categories.
const analyticBaseApiCost = 18e9 / 1e6 * 0.03 + 6e9 / 1e6 * 1.5 + 6e9 / 1e6 * 6;
near(baseline.firstMonth.fullApiCost, analyticBaseApiCost, 'Month-one API token bill');
near(baseline.NPV, -defaults.capex + (analyticBaseApiCost * defaults.availability - 24_000) * annuityFactor + residualPV, 'Base NPV');
near(thresholds.monthlyOperatingBreakEven.demandMultiplier, defaults.monthlyOpex / uncongestedAvoidedApiCost, 'Operating threshold');
const peakCapacity = baseline.capacity.buckets[0];
const cappedPeakAvoidedApi = defaults.availability * defaults.headroom * defaults.peakHours * 3600 * (defaults.inputTokensPerSecond * 0.53 + defaults.outputTokensPerSecond * 8) / 1e6;
const offPeakAvoidedApiPerDemandMultiplier = defaults.availability * (12e9 * 0.265 + 3e9 * 4) / 1e6;
const requiredAvoidedCost = defaults.monthlyOpex + (defaults.capex - residualPV) / annuityFactor;
const analyticNpvThreshold = (requiredAvoidedCost - cappedPeakAvoidedApi) / offPeakAvoidedApiPerDemandMultiplier;
const analyticOperatingPaybackThreshold = (defaults.monthlyOpex + defaults.capex / defaults.horizon - cappedPeakAvoidedApi) / offPeakAvoidedApiPerDemandMultiplier;
const analyticDiscountedOperatingPaybackThreshold = (defaults.monthlyOpex + defaults.capex / annuityFactor - cappedPeakAvoidedApi) / offPeakAvoidedApiPerDemandMultiplier;
assert.ok(analyticNpvThreshold > peakCapacity.noCongestionDemandMultiplier && analyticNpvThreshold < saturationDemandMultiplier);
near(thresholds.npvBreakEvenIncludingTerminal.demandMultiplier, analyticNpvThreshold, 'Capacity-aware NPV threshold');
near(thresholds.operatingPaybackWithinHorizon.demandMultiplier, analyticOperatingPaybackThreshold, 'Operating payback threshold');
near(thresholds.discountedOperatingPaybackWithinHorizon.demandMultiplier, analyticDiscountedOperatingPaybackThreshold, 'Discounted operating payback threshold');

const noResale = evaluate({ terminalResidual: 0 });
const resaleFraction = 0.15;
const maxInitialInvestmentWithProportionalResale = baseline.presentValueOperatingSavings / (1 - resaleFraction * terminalDiscountFactor);
const atJointResaleLimit = evaluate({ capex: maxInitialInvestmentWithProportionalResale, terminalResidual: resaleFraction * maxInitialInvestmentWithProportionalResale });
near(baseline.maxInitialInvestment, 20_629.2 * annuityFactor + residualPV, 'Fixed-resale investment limit');
near(noResale.maxInitialInvestment, 20_629.2 * annuityFactor, 'Zero-resale investment limit');
near(evaluate({ capex: baseline.maxInitialInvestment }).NPV, 0, 'NPV at fixed-resale limit');
near(atJointResaleLimit.NPV, 0, 'NPV at proportional-resale limit');
near(baseline.idealServiceCeiling, 45_540 * annuityFactor + residualPV, 'Ideal service ceiling');
const doubleCapacity = evaluate({ capacityMultiplier: 2 });
const mid = evaluate(midDemand);
const midMoreCapacity = evaluate({ ...midDemand, capacityMultiplier: 1.25 });
near(doubleCapacity.NPV, baseline.NPV, 'Uncongested demand cannot value extra throughput');
near(midMoreCapacity.NPV - mid.NPV, cappedPeakAvoidedApi * 0.25 * annuityFactor, 'Additional simultaneous capacity value at congested peak');
const saturated = evaluate(demandAt(saturationDemandMultiplier));
const moreDemand = evaluate(demandAt(saturationDemandMultiplier * 2));
const analyticMaxAvoidedApi = defaults.availability * defaults.headroom * 3600 / 1e6 * (
  defaults.peakHours * (defaults.inputTokensPerSecond * 0.53 + defaults.outputTokensPerSecond * 8)
  + defaults.offHours * (defaults.inputTokensPerSecond * 0.265 + defaults.outputTokensPerSecond * 4));
const analyticMaxMonthlySavings = analyticMaxAvoidedApi - defaults.monthlyOpex;
near(saturated.firstMonth.avoidedApiCost, analyticMaxAvoidedApi, 'Saturated avoided-API envelope');
near(saturated.NPV, -defaults.capex + analyticMaxMonthlySavings * annuityFactor + residualPV, 'Saturated NPV envelope');
near(moreDemand.NPV, saturated.NPV, 'Same-cache saving plateau');
const high = evaluate(highDemand);
const coldFallback = evaluate({ ...highDemand, fallbackCacheHit: 0 });
const analyticCachePenalty = high.firstMonth.buckets.reduce((sum, b) => sum + b.overflowInputTokens * 0.75 * (defaults.prices[b.id].miss - defaults.prices[b.id].hit) / 1e6, 0);
near(coldFallback.firstMonth.cachePenalty, analyticCachePenalty, 'Independent fallback input-cache penalty');
near(coldFallback.NPV, high.NPV - analyticCachePenalty * annuityFactor, 'Fallback-cache cash-flow effect');
const firstPrinciples = {
  capitalLimits: {
    maxInitialInvestment: baseline.maxInitialInvestment, maxInitialInvestmentWithoutResale: baseline.maxInitialInvestmentWithoutResale,
    idealServiceCeiling: baseline.idealServiceCeiling, fixedTerminalResidual: defaults.terminalResidual, resaleFraction,
    maxInitialInvestmentWithProportionalResale, jointTerminalResidual: resaleFraction * maxInitialInvestmentWithProportionalResale,
    zeroOpexNPV: evaluate({ monthlyOpex: 0 }).NPV,
  },
  throughput: {
    referenceInputTokensPerSecond: defaults.inputTokensPerSecond, referenceOutputTokensPerSecond: defaults.outputTokensPerSecond,
    baseCapacity1To2MonthlySavingsGain: doubleCapacity.firstMonth.operatingSavings - baseline.firstMonth.operatingSavings,
    baseCapacity1To2NPVGain: doubleCapacity.NPV - baseline.NPV,
    baseMinimumUncongestedCapacityMultiplier: 1 / baseline.capacity.maximumNoCongestionDemandMultiplier,
    midCapacity1To1_25MonthlySavingsGain: midMoreCapacity.firstMonth.operatingSavings - mid.firstMonth.operatingSavings,
    midCapacity1To1_25NPVGain: midMoreCapacity.NPV - mid.NPV,
    midMinimumUncongestedCapacityMultiplier: 1 / mid.capacity.maximumNoCongestionDemandMultiplier,
  },
  saturationEnvelope: {
    assumptions: 'Fixed simultaneous token capacities, 4:1 input/output mix, 50/50 bucket splits, default prices and costs, equal baseline/fallback cache hit rates, sufficient external API capacity.',
    ...thresholdAt(saturationDemandMultiplier), localInputTokens: saturated.firstMonth.localInputTokens, localOutputTokens: saturated.firstMonth.localOutputTokens,
    avoidedApiCost: saturated.firstMonth.avoidedApiCost, operatingSavings: saturated.firstMonth.operatingSavings,
    maxInitialInvestment: saturated.maxInitialInvestment, NPV: saturated.NPV, capitalROI: saturated.capitalROI,
    simplePaybackMonth: saturated.simplePaybackMonth, discountedPaybackMonth: saturated.discountedPaybackMonth,
  },
  fallbackCacheLoss: {
    ...highDemand, apiCacheHit: defaults.apiCacheHit, fallbackCacheHit: 0,
    monthlyPenalty: coldFallback.firstMonth.cachePenalty, presentValuePenalty: analyticCachePenalty * annuityFactor, NPV: coldFallback.NPV,
  },
};

// R includes all recurring rental-service costs, replaces ownership OpEx, and excludes fallback.
const rentalSetup = 100_000;
function rentalCase(demand) {
  const r = evaluate(demand);
  const fallbackPV = r.months.reduce((sum, m) => sum + m.overflowApiCost * (m.discountFactor ?? 0), 0);
  const avoidedApiPV = r.presentValueFullApi - fallbackPV;
  const maxMonthlyRentalRate = (avoidedApiPV - rentalSetup) / annuityFactor;
  near(evaluate({ ...demand, capex: rentalSetup, monthlyOpex: maxMonthlyRentalRate, terminalResidual: 0 }).NPV, 0, 'Rental NPV at rate limit');
  near(maxMonthlyRentalRate, r.firstMonth.avoidedApiCost - rentalSetup / annuityFactor, 'Constant-flow rental limit');
  return {
    ...demand, allApiTco: r.fullApiTco, ownershipTco: r.hybridTco,
    allApiPresentValueCost: r.presentValueFullApi, ownershipPresentValueCost: r.presentValueFullApi - r.NPV,
    ownershipNPV: r.NPV, monthlyFallbackApiCost: r.firstMonth.overflowApiCost, monthlyAvoidedApiCost: r.firstMonth.avoidedApiCost,
    rentalTcoIntercept: rentalSetup + r.totalOverflowApiCost, rentalTcoRateCoefficient: defaults.horizon,
    rentalPresentValueCostIntercept: rentalSetup + fallbackPV, rentalPresentValueCostRateCoefficient: annuityFactor,
    rentalNPVIntercept: avoidedApiPV - rentalSetup, maxMonthlyRentalRate, maxMonthlyRentalRateZeroSetup: avoidedApiPV / annuityFactor,
  };
}
const scaledPrices = factor => Object.fromEntries(Object.entries(defaults.prices).map(([bucket, prices]) => [bucket, Object.fromEntries(Object.entries(prices).map(([kind, price]) => [kind, price * factor]))]));
const priceMultiplierNPV = factor => evaluate({ ...highDemand, prices: scaledPrices(factor) }).NPV;
const priceMultiplierBoundary = bisectIncreasing(priceMultiplierNPV, 0, 1);
const immediatePriceReduction = 1 - priceMultiplierBoundary;
const capacityNPV = capacityMultiplier => evaluate({ ...highDemand, capacityMultiplier }).NPV;
const minimumCapacityMultiplier = bisectIncreasing(capacityNPV, 0, 1);
const fallbackCacheNPV = fallbackCacheHit => evaluate({ ...highDemand, fallbackCacheHit }).NPV;
const minimumFallbackCacheHit = bisectIncreasing(fallbackCacheNPV, 0, defaults.apiCacheHit);
const analyticPriceMultiplierBoundary = requiredAvoidedCost / high.firstMonth.avoidedApiCost;
const highOffPeakAvoidedApi = defaults.availability * (32e9 * 0.265 + 8e9 * 4) / 1e6;
const analyticMinimumCapacityMultiplier = (requiredAvoidedCost - highOffPeakAvoidedApi) / cappedPeakAvoidedApi;
const cachePenaltyPerHitFraction = high.firstMonth.buckets.reduce((sum, b) => sum + b.overflowInputTokens * (defaults.prices[b.id].miss - defaults.prices[b.id].hit) / 1e6, 0);
const analyticMinimumFallbackCacheHit = defaults.apiCacheHit - high.NPV / (cachePenaltyPerHitFraction * annuityFactor);
near(priceMultiplierBoundary, analyticPriceMultiplierBoundary, 'Immediate price boundary');
near(minimumCapacityMultiplier, analyticMinimumCapacityMultiplier, 'Simultaneous capacity boundary');
near(minimumFallbackCacheHit, analyticMinimumFallbackCacheHit, 'Fallback cache boundary');
const atCapacityBoundary = evaluate({ ...highDemand, capacityMultiplier: minimumCapacityMultiplier });
assert.ok(atCapacityBoundary.firstMonth.buckets[0].onlineFraction < 1 && atCapacityBoundary.firstMonth.buckets[1].onlineFraction === 1);
const probeStep = 0.000001;
function boundaryProbe(value, fn, increasing = true) {
  const result = { value, NPVatBoundary: fn(value), probeStep, belowValue: value - probeStep, NPVbelow: fn(value - probeStep), aboveValue: value + probeStep, NPVabove: fn(value + probeStep) };
  near(result.NPVatBoundary, 0, 'NPV at boundary');
  assert.ok(increasing ? result.NPVbelow < 0 && result.NPVabove > 0 : result.NPVbelow > 0 && result.NPVabove < 0, 'NPV changes sign across boundary');
  return result;
}
function stressCase(overrides) {
  const effectiveOverrides = { ...highDemand, ...overrides };
  const r = evaluate(effectiveOverrides);
  return { overrides: effectiveOverrides, NPV: r.NPV, fullApiTco: r.fullApiTco, hybridTco: r.hybridTco,
    firstMonthOperatingSavings: r.firstMonth.operatingSavings, lastMonthOperatingSavings: r.months.at(-1).operatingSavings,
    capitalROI: r.capitalROI, simplePaybackMonth: r.simplePaybackMonth, discountedPaybackMonth: r.discountedPaybackMonth };
}
const apiDiscountedMultiplierSum = Array.from({ length: defaults.horizon }, (_, i) => 0.7 ** (i / 12) * 1.1 ** (-(i + 1) / 12)).reduce((sum, x) => sum + x, 0);
const reducedDemand = { monthlyInputTokens: 51.2e9, monthlyOutputTokens: 12.8e9 };
const apiAnnualMinus30 = stressCase({ apiAnnualChange: -0.3 });
const monthlyDemandMinus20 = stressCase(reducedDemand);
const combined = stressCase({ ...reducedDemand, apiAnnualChange: -0.3 });
const reducedDemandAvoidedApi = cappedPeakAvoidedApi + defaults.availability * (25.6e9 * 0.265 + 6.4e9 * 4) / 1e6;
near(apiAnnualMinus30.NPV, -defaults.capex + high.firstMonth.avoidedApiCost * apiDiscountedMultiplierSum - defaults.monthlyOpex * annuityFactor + residualPV, 'High-demand annual API reduction');
near(monthlyDemandMinus20.NPV, -defaults.capex + (reducedDemandAvoidedApi - defaults.monthlyOpex) * annuityFactor + residualPV, 'High demand reduced once by 20%');
near(combined.NPV, -defaults.capex + reducedDemandAvoidedApi * apiDiscountedMultiplierSum - defaults.monthlyOpex * annuityFactor + residualPV, 'Combined stress');
const decisionMemo = {
  rental: {
    currency: 'RMB', horizonMonths: defaults.horizon, annualDiscount: defaults.annualDiscount, oneTimeSetup: rentalSetup, terminalResidual: 0,
    monthlyRateDefinition: 'R is complete recurring rental-service cost: node, energy, network, storage, maintenance, software and operations; excludes separately billed fallback. Ownership OpEx is replaced, never added again.',
    serviceAssumptions: 'Same simultaneously attainable token rates, availability, token workload, quality, routing and fallback as ownership. Conditional on service validation; not a quote.',
    base: rentalCase(demandAt(1)), high: rentalCase(highDemand),
    monthlyRateEqualToOwnershipPresentValue: defaults.monthlyOpex + (defaults.capex - rentalSetup - residualPV) / annuityFactor,
  },
  fragilityHigh: {
    ...highDemand, baselineNPV: high.NPV,
    assumptions: 'One-at-a-time changes; all other default parameters fixed. Zero-NPV boundaries are not statistically justified safety margins.',
    immediateApiPriceReduction: { ...boundaryProbe(immediatePriceReduction, reduction => priceMultiplierNPV(1 - reduction), false), percent: immediatePriceReduction * 100, remainingPriceMultiplier: priceMultiplierBoundary, definition: 'Every hit, miss and output tariff in both buckets falls immediately by the same fraction and stays flat through the horizon; both baseline and fallback bills change.' },
    minimumCapacityMultiplier: { ...boundaryProbe(minimumCapacityMultiplier, capacityNPV), inputTokensPerSecond: defaults.inputTokensPerSecond * minimumCapacityMultiplier, outputTokensPerSecond: defaults.outputTokensPerSecond * minimumCapacityMultiplier, definition: 'Scales both jointly attainable rates before headroom and availability; not independent peak capacities.' },
    minimumFallbackCacheHit: { ...boundaryProbe(minimumFallbackCacheHit, fallbackCacheNPV), percent: minimumFallbackCacheHit * 100, definition: 'Input-token-weighted fallback cache hit fraction; all-API baseline stays at 75%.' },
  },
  stressesHigh: { definitions: 'API prices compound at 0.7^((month-1)/12). Demand reduction is one change to 51.2B input / 12.8B output per month, not repeated monthly decline.', apiAnnualMinus30, monthlyDemandMinus20, combined },
  independentChecks: { analyticPriceMultiplierBoundary, analyticMinimumCapacityMultiplier, analyticMinimumFallbackCacheHit, apiDiscountedMultiplierSum, passed: true },
};
const document = {
  schemaVersion: 5,
  engine: {
    financeModule: '../llm-economics/finance.js',
    tokenModule: '../llm-economics/token-model.js',
    profileModule: '../llm-economics/presets.js',
    scope: 'Generic month-end investment arithmetic; this explicit profile applies a token-based API-cost-substitution scenario.',
  },
  profile: examplePreset.metadata,
  assumptions: {
    currency: 'RMB', tokenPriceUnit: 'RMB per million tokens', demandUnit: 'Input and output tokens per month, separately specified', cashFlow: 'Pre-tax, month-end; no depreciation expense',
    status: 'Illustrative budgets and capacities, not a quote or hardware benchmark',
    month: '30 days; peak 154h = 22 illustrative non-holiday workdays × 7h; off-peak 566h',
    capacity: 'Input/output rates must be simultaneously achieved for the specified workload and SLO. Proportional bucket routing uses the tighter token constraint and cannot predict real scheduling or bursts. CapacityMultiplier scales both rates together.',
    availability: 'Uniform arrivals independent of downtime: f = a × min(1, inputCapacity/inputDemand, outputCapacity/outputDemand). Ignore zero-demand dimensions; both zero process nothing. Downtime traffic falls back immediately; headroom is applied once.',
    unconstrainedThresholds: 'Hold 4:1 input/output mix and bucket shares fixed; ignore congestion but retain downtime fallback. Demand multipliers reference 24B input / 6B output monthly.',
    apiCache: 'Input-token-weighted API hit fraction affects billing, never local input-capacity demand. It is not local prefix caching.',
    fallbackCache: 'Null follows apiCacheHit; an explicit value prices fallback independently. cachePenalty is the signed extra fallback bill relative to baseline caching.',
    variableCost: 'Per-million local input/output costs are additional to monthlyOpex; do not duplicate electricity or labor already included. opexAnnualGrowth scales fixed costs and these per-token rates. Routing stays local-first.',
    output: 'Includes all billed thinking output once; no image workload assumed',
    capitalAllocation: '(capex − residual)/horizon and capex/horizon are comparison references, not extra cash expenses',
    residual: 'Final-month expected sale proceeds; included in terminal cash/NPV, never operating payback',
    capexSensitivity: 'Optional initial-investment reduction keeps absolute residual fixed',
    investmentLimits: 'maxInitialInvestment holds residual fixed; proportional resale is solved jointly. idealServiceCeiling assumes zero service costs/fallback and is not attainable revenue.',
    scope: 'No financing, taxes, replacement CapEx, price elasticity or service-quality change modeled. Horizon and independent input/output/cost growth are adjustable.',
  },
  defaults, thresholds, firstPrinciples, decisionMemo,
  independentChecks: { annuityFactor, residualPV, analyticBaseApiCost, analyticNpvThreshold, analyticOperatingPaybackThreshold, analyticDiscountedOperatingPaybackThreshold, analyticMaxAvoidedApi, analyticMaxMonthlySavings, analyticCachePenalty, proportionalResaleLimitNPV: atJointResaleLimit.NPV, passed: true },
  scenarios: definitions.map(definition => ({ ...definition, result: evaluate(definition.overrides) })),
};
// Exercise the standalone bill helper as well as evaluation, without maintaining a second pricing engine.
near(billTokens(12e9, 3e9).total + billTokens(12e9, 3e9, {}, 'off').total, analyticBaseApiCost, 'Public token bill helper');
const json = `${JSON.stringify(document, null, 2)}\n`;
if (process.argv.includes('--stdout')) process.stdout.write(json);
else if (process.argv.includes('--check')) {
  assert.equal(await readFile(destination, 'utf8'), json, 'scenarios.json is stale; run node scripts/calculate-deepseek-roi.mjs');
  console.log('DeepSeek token ROI scenarios and independent arithmetic checks passed.');
} else {
  await writeFile(destination, json);
  console.log(`Wrote ${destination}`);
  for (const { id, result } of document.scenarios) console.log(`${id}: API TCO ¥${result.fullApiTco.toFixed(2)}, hybrid TCO ¥${result.hybridTco.toFixed(2)}, NPV ¥${result.NPV.toFixed(2)}, operating payback ${result.simplePaybackMonth ?? 'none'} months`);
}
