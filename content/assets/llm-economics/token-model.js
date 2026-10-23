/**
 * Currency-neutral token service-cost comparison against a full-API baseline.
 * Every parameter and tariff must be explicitly supplied in one currency; no
 * vendor preset or metadata is read by this module. Prices/cost rates are per
 * million tokens. Input/output rates must be jointly attainable under the same
 * workload and SLO, not separate peak benchmarks. Bucket streams share one routing
 * fraction, making this a screening model rather than a scheduling prediction.
 * API cache hits change bills, never local token demand. Uniform arrivals are
 * independent of downtime; traffic during downtime immediately falls back.
 */
import { evaluateCashFlows } from './finance.js';

const parameterKeys = [
  'monthlyInputTokens', 'monthlyOutputTokens', 'inputTokensPerSecond', 'outputTokensPerSecond',
  'capacityMultiplier', 'inputPeakShare', 'outputPeakShare', 'inputAnnualGrowth', 'outputAnnualGrowth',
  'opexAnnualGrowth', 'perMillionLocalInputCost', 'perMillionLocalOutputCost', 'apiCacheHit',
  'fallbackCacheHit', 'peakHours', 'offHours', 'availability', 'headroom', 'capex', 'monthlyOpex',
  'horizon', 'annualDiscount', 'terminalResidual', 'apiAnnualChange', 'prices',
];

function validateParameters(parameters) {
  if (parameters === null || typeof parameters !== 'object' || Array.isArray(parameters)) throw new RangeError('A complete explicit token profile is required');
  for (const key of Object.keys(parameters)) {
    if (!parameterKeys.includes(key)) throw new RangeError(`Unknown token-model parameter: ${key}`);
  }
  for (const key of parameterKeys) {
    if (!Object.hasOwn(parameters, key)) throw new RangeError(`Missing token-model parameter: ${key}`);
  }
  for (const bucket of ['peak', 'off']) {
    for (const type of ['hit', 'miss', 'output']) {
      if (!Object.hasOwn(parameters.prices?.[bucket] ?? {}, type)) throw new RangeError(`Missing prices.${bucket}.${type}`);
    }
  }
  const p = {
    ...parameters,
    prices: { peak: { ...parameters.prices.peak }, off: { ...parameters.prices.off } },
  };
  for (const key of [
    'monthlyInputTokens', 'monthlyOutputTokens', 'inputTokensPerSecond',
    'outputTokensPerSecond', 'capacityMultiplier', 'peakHours', 'offHours',
    'capex', 'monthlyOpex', 'perMillionLocalInputCost', 'perMillionLocalOutputCost', 'terminalResidual',
  ]) {
    if (!Number.isFinite(p[key]) || p[key] < 0) throw new RangeError(`${key} must be finite and nonnegative`);
  }
  for (const key of ['apiCacheHit', 'inputPeakShare', 'outputPeakShare', 'availability', 'headroom']) {
    if (!Number.isFinite(p[key]) || p[key] < 0 || p[key] > 1) throw new RangeError(`${key} must be between 0 and 1`);
  }
  if (p.fallbackCacheHit !== null && (!Number.isFinite(p.fallbackCacheHit) || p.fallbackCacheHit < 0 || p.fallbackCacheHit > 1)) {
    throw new RangeError('fallbackCacheHit must be null or between 0 and 1');
  }
  for (const key of ['annualDiscount', 'apiAnnualChange', 'inputAnnualGrowth', 'outputAnnualGrowth', 'opexAnnualGrowth']) {
    if (!Number.isFinite(p[key]) || p[key] <= -1) throw new RangeError(`${key} must be finite and greater than -1`);
  }
  if (!Number.isInteger(p.horizon) || p.horizon < 1) throw new RangeError('horizon must be a positive integer');
  if (p.peakHours + p.offHours <= 0) throw new RangeError('The month must have positive hours');
  for (const bucket of ['peak', 'off']) {
    for (const type of ['hit', 'miss', 'output']) {
      if (!Number.isFinite(p.prices[bucket][type]) || p.prices[bucket][type] < 0) throw new RangeError(`prices.${bucket}.${type} must be finite and nonnegative`);
    }
  }
  return p;
}

function tokenBill(inputCount, outputCount, p, bucket, month, route) {
  const hit = route === 'fallback' ? (p.fallbackCacheHit ?? p.apiCacheHit) : p.apiCacheHit;
  const factor = (1 + p.apiAnnualChange) ** ((month - 1) / 12) / 1_000_000;
  const counts = { hit: inputCount * hit, miss: inputCount * (1 - hit), output: outputCount };
  const result = Object.fromEntries(Object.entries(counts).map(([type, tokens]) => [type, { tokens, amount: tokens * p.prices[bucket][type] * factor }]));
  return { ...result, total: result.hit.amount + result.miss.amount + result.output.amount };
}

/** Three billed token categories, each with token count and amount, plus total. */
export function billTokens(inputTokenCount, outputTokenCount, parameters, bucket = 'peak', month = 1, route = 'baseline') {
  for (const count of [inputTokenCount, outputTokenCount]) {
    if (!Number.isFinite(count) || count < 0) throw new RangeError('Token counts must be finite and nonnegative');
  }
  if (!['peak', 'off'].includes(bucket)) throw new RangeError('bucket must be peak or off');
  if (!Number.isInteger(month) || month < 1) throw new RangeError('month must be a positive integer');
  if (!['baseline', 'fallback'].includes(route)) throw new RangeError('route must be baseline or fallback');
  return tokenBill(inputTokenCount, outputTokenCount, validateParameters(parameters), bucket, month, route);
}

function sumBills(bills) {
  const result = { hit: { tokens: 0, amount: 0 }, miss: { tokens: 0, amount: 0 }, output: { tokens: 0, amount: 0 }, total: 0 };
  for (const bill of bills) {
    for (const type of ['hit', 'miss', 'output']) {
      result[type].tokens += bill[type].tokens;
      result[type].amount += bill[type].amount;
    }
    result.total += bill.total;
  }
  return result;
}

/**
 * Month zero pays CapEx; later cash flows occur at month end. Residual is received
 * only in the final month and never creates an operating-payback crossing.
 * Fixed and additional local per-token operating costs grow at opexAnnualGrowth.
 * Routing stays local-first; no dispatch optimization, financing, taxes, hardware
 * replacement, or improvement in task quality is silently assumed.
 */
export function evaluateTokenScenario(parameters) {
  const p = validateParameters(parameters);
  const inputRate = p.inputTokensPerSecond * p.capacityMultiplier;
  const outputRate = p.outputTokensPerSecond * p.capacityMultiplier;
  const definitions = [
    { id: 'peak', hours: p.peakHours, inputShare: p.inputPeakShare, outputShare: p.outputPeakShare },
    { id: 'off', hours: p.offHours, inputShare: 1 - p.inputPeakShare, outputShare: 1 - p.outputPeakShare },
  ];
  const capacityBuckets = definitions.map(b => {
    const onlineInputTokenCapacity = b.hours * 3600 * p.headroom * inputRate;
    const onlineOutputTokenCapacity = b.hours * 3600 * p.headroom * outputRate;
    const firstInput = p.monthlyInputTokens * b.inputShare;
    const firstOutput = p.monthlyOutputTokens * b.outputShare;
    const constraints = [];
    if (firstInput > 0) constraints.push(onlineInputTokenCapacity / firstInput);
    if (firstOutput > 0) constraints.push(onlineOutputTokenCapacity / firstOutput);
    return {
      ...b, onlineInputTokenCapacity, onlineOutputTokenCapacity,
      effectiveInputTokenCapacity: p.availability * onlineInputTokenCapacity,
      effectiveOutputTokenCapacity: p.availability * onlineOutputTokenCapacity,
      noCongestionDemandMultiplier: constraints.length ? Math.min(...constraints) : null,
    };
  });
  const activeLimits = capacityBuckets.map(b => b.noCongestionDemandMultiplier).filter(x => x !== null);
  const maximumNoCongestionDemandMultiplier = activeLimits.length ? Math.min(...activeLimits) : null;
  const capacity = {
    capacityMultiplier: p.capacityMultiplier,
    assumedInputTokensPerSecond: inputRate, assumedOutputTokensPerSecond: outputRate,
    onlineInputTokensPerSecond: inputRate * p.headroom,
    onlineOutputTokensPerSecond: outputRate * p.headroom,
    effectiveInputTokensPerSecond: inputRate * p.headroom * p.availability,
    effectiveOutputTokensPerSecond: outputRate * p.headroom * p.availability,
    totalMonthlyInputTokenCapacity: capacityBuckets.reduce((sum, b) => sum + b.effectiveInputTokenCapacity, 0),
    totalMonthlyOutputTokenCapacity: capacityBuckets.reduce((sum, b) => sum + b.effectiveOutputTokenCapacity, 0),
    maximumNoCongestionDemandMultiplier,
    maximumAllLocalDemandMultiplier: p.availability < 1 ? 0 : maximumNoCongestionDemandMultiplier,
    buckets: capacityBuckets,
  };
  const tokenFields = [
    'demandInputTokens', 'demandOutputTokens', 'localInputTokens', 'localOutputTokens',
    'overflowInputTokens', 'overflowOutputTokens', 'downtimeOverflowInputTokens',
    'downtimeOverflowOutputTokens', 'overloadOverflowInputTokens', 'overloadOverflowOutputTokens',
  ];
  const zeroTokens = Object.fromEntries(tokenFields.map(key => [key, 0]));
  const serviceMonths = [{
    month: 0, ...zeroTokens, localInputShare: null, localOutputShare: null,
    fullApiCost: 0, overflowApiCost: 0, fixedOpex: 0, variableLocalCost: 0,
    variableLocalInputCost: 0, variableLocalOutputCost: 0, opex: 0, hybridOperatingCost: 0,
    localDisplacedApiCost: 0, cachePenalty: 0, avoidedApiCost: 0, operatingSavings: 0, residual: 0,
    baselineBillBreakdown: sumBills([]), fallbackBillBreakdown: sumBills([]),
  }];
  const totals = { ...zeroTokens, fullApiCost: 0, overflowApiCost: 0, fixedOpex: 0, variableLocalCost: 0, localDisplacedApiCost: 0, cachePenalty: 0 };
  let presentValueFullApi = 0;
  for (let month = 1; month <= p.horizon; month++) {
    const demandInputTokens = p.monthlyInputTokens * (1 + p.inputAnnualGrowth) ** ((month - 1) / 12);
    const demandOutputTokens = p.monthlyOutputTokens * (1 + p.outputAnnualGrowth) ** ((month - 1) / 12);
    const opexFactor = (1 + p.opexAnnualGrowth) ** ((month - 1) / 12);
    const fixedOpex = p.monthlyOpex * opexFactor;
    if (![demandInputTokens, demandOutputTokens, opexFactor, fixedOpex].every(Number.isFinite)) throw new RangeError('Projected demand or operating cost is not finite');
    const buckets = capacityBuckets.map(b => {
      const input = demandInputTokens * b.inputShare;
      const output = demandOutputTokens * b.outputShare;
      const constraints = [1];
      if (input > 0) constraints.push(b.onlineInputTokenCapacity / input);
      if (output > 0) constraints.push(b.onlineOutputTokenCapacity / output);
      const onlineFraction = input > 0 || output > 0 ? Math.min(...constraints) : 0;
      const localFraction = p.availability * onlineFraction;
      const localInputTokens = input * localFraction;
      const localOutputTokens = output * localFraction;
      const overflowInputTokens = input - localInputTokens;
      const overflowOutputTokens = output - localOutputTokens;
      const baselineBillBreakdown = tokenBill(input, output, p, b.id, month, 'baseline');
      const fallbackBillBreakdown = tokenBill(overflowInputTokens, overflowOutputTokens, p, b.id, month, 'fallback');
      const unchangedFallbackBill = tokenBill(overflowInputTokens, overflowOutputTokens, p, b.id, month, 'baseline');
      const localDisplacedApiCost = tokenBill(localInputTokens, localOutputTokens, p, b.id, month, 'baseline').total;
      const cachePenalty = fallbackBillBreakdown.total - unchangedFallbackBill.total;
      const variableLocalInputCost = localInputTokens / 1_000_000 * p.perMillionLocalInputCost * opexFactor;
      const variableLocalOutputCost = localOutputTokens / 1_000_000 * p.perMillionLocalOutputCost * opexFactor;
      return {
        ...b, onlineFraction, localFraction, demandInputTokens: input, demandOutputTokens: output,
        localInputTokens, localOutputTokens, overflowInputTokens, overflowOutputTokens,
        downtimeOverflowInputTokens: (1 - p.availability) * input,
        downtimeOverflowOutputTokens: (1 - p.availability) * output,
        overloadOverflowInputTokens: p.availability * (1 - onlineFraction) * input,
        overloadOverflowOutputTokens: p.availability * (1 - onlineFraction) * output,
        localInputShare: input > 0 ? localFraction : null,
        localOutputShare: output > 0 ? localFraction : null,
        baselineBillBreakdown, fallbackBillBreakdown,
        fullApiCost: baselineBillBreakdown.total, overflowApiCost: fallbackBillBreakdown.total,
        localDisplacedApiCost, cachePenalty,
        avoidedApiCost: baselineBillBreakdown.total - fallbackBillBreakdown.total,
        variableLocalInputCost, variableLocalOutputCost,
        variableLocalCost: variableLocalInputCost + variableLocalOutputCost,
      };
    });
    const sum = key => buckets.reduce((total, b) => total + b[key], 0);
    const tokens = Object.fromEntries(tokenFields.map(key => [key, sum(key)]));
    const fullApiCost = sum('fullApiCost');
    const overflowApiCost = sum('overflowApiCost');
    const variableLocalCost = sum('variableLocalCost');
    const opex = fixedOpex + variableLocalCost;
    const hybridOperatingCost = opex + overflowApiCost;
    const operatingSavings = fullApiCost - hybridOperatingCost;
    const discountFactor = (1 + p.annualDiscount) ** (-month / 12);
    const row = {
      month, ...tokens,
      localInputShare: tokens.demandInputTokens > 0 ? tokens.localInputTokens / tokens.demandInputTokens : null,
      localOutputShare: tokens.demandOutputTokens > 0 ? tokens.localOutputTokens / tokens.demandOutputTokens : null,
      fullApiCost, overflowApiCost, fixedOpex, variableLocalCost,
      variableLocalInputCost: sum('variableLocalInputCost'), variableLocalOutputCost: sum('variableLocalOutputCost'),
      opex, hybridOperatingCost, localDisplacedApiCost: sum('localDisplacedApiCost'), cachePenalty: sum('cachePenalty'),
      avoidedApiCost: fullApiCost - overflowApiCost, operatingSavings,
      baselineBillBreakdown: sumBills(buckets.map(b => b.baselineBillBreakdown)),
      fallbackBillBreakdown: sumBills(buckets.map(b => b.fallbackBillBreakdown)),
    };
    serviceMonths.push(month === 1 ? { ...row, buckets } : row);
    for (const key of Object.keys(totals)) totals[key] += row[key];
    presentValueFullApi += fullApiCost * discountFactor;
  }
  const totalOpex = totals.fixedOpex + totals.variableLocalCost;
  const hybridTco = p.capex + totalOpex + totals.overflowApiCost - p.terminalResidual;
  const finance = evaluateCashFlows({
    initialInvestment: p.capex,
    monthlyNetCashFlows: serviceMonths.slice(1).map(row => row.operatingSavings),
    terminalValue: p.terminalResidual,
    annualDiscount: p.annualDiscount,
  });
  const months = serviceMonths.map((service, month) => {
    const { buckets, ...row } = service;
    const cash = finance.months[month];
    return {
      ...row,
      ...(month > 0 ? { discountFactor: cash.discountFactor } : {}),
      residual: cash.terminalValue, cashFlow: cash.cashFlow, discountedCashFlow: cash.discountedCashFlow,
      cumulativeOperatingNet: cash.cumulativePreTerminalNet,
      cumulativeDiscountedOperatingNet: cash.cumulativeDiscountedPreTerminalNet,
      cumulativeNet: cash.cumulativeNet, cumulativeDiscountedNet: cash.cumulativeDiscountedNet,
    };
  });
  const firstMonth = {
    ...months[1], buckets: serviceMonths[1].buckets,
    capitalAllocatedCost: (p.capex - p.terminalResidual) / p.horizon,
    grossCapitalAllocatedCost: p.capex / p.horizon,
    allocatedHybridCost: months[1].hybridOperatingCost + (p.capex - p.terminalResidual) / p.horizon,
    grossAllocatedHybridCost: months[1].hybridOperatingCost + p.capex / p.horizon,
  };
  const presentValueOperatingSavings = finance.presentValueMonthlyNetCashFlows;
  const presentValueResidual = finance.presentValueTerminalValue;
  return {
    parameters: p, capacity, monthlyDiscount: finance.monthlyDiscount,
    firstMonth, months, fullApiTco: totals.fullApiCost, hybridTco,
    totalOpex, totalFixedOpex: totals.fixedOpex, totalVariableLocalCost: totals.variableLocalCost,
    totalOverflowApiCost: totals.overflowApiCost, totalLocalDisplacedApiCost: totals.localDisplacedApiCost,
    totalCachePenalty: totals.cachePenalty,
    totalDemandInputTokens: totals.demandInputTokens, totalDemandOutputTokens: totals.demandOutputTokens,
    totalLocalInputTokens: totals.localInputTokens, totalLocalOutputTokens: totals.localOutputTokens,
    totalOverflowInputTokens: totals.overflowInputTokens, totalOverflowOutputTokens: totals.overflowOutputTokens,
    presentValueFullApi, presentValueOperatingSavings, presentValueResidual,
    maxInitialInvestment: presentValueOperatingSavings + presentValueResidual,
    maxInitialInvestmentWithoutResale: presentValueOperatingSavings,
    idealServiceCeiling: presentValueFullApi + presentValueResidual,
    operatingNetCash: finance.preTerminalNetCash, terminalNetCash: finance.terminalNetCash,
    NPV: finance.NPV, operatingNPV: finance.preTerminalNPV,
    capitalROI: finance.capitalROI,
    simplePaybackMonth: finance.simplePaybackMonth, discountedPaybackMonth: finance.discountedPaybackMonth,
    simplePaybackReversalMonths: finance.simplePaybackReversalMonths, discountedPaybackReversalMonths: finance.discountedPaybackReversalMonths,
    sustainedSimplePaybackMonth: finance.sustainedSimplePaybackMonth, sustainedDiscountedPaybackMonth: finance.sustainedDiscountedPaybackMonth,
  };
}
