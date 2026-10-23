/** Explicit profiles; none is silently selected by the generic model.
 * Metadata identifies an example or a user-supplied scenario, never changes math.
 * Currency is a unit label, not an exchange-rate conversion. Every amount in a
 * profile must already be expressed in that one currency.
 */
const exampleParameters = Object.freeze({
  monthlyInputTokens: 24_000_000_000,
  monthlyOutputTokens: 6_000_000_000, // Includes all billed thinking output once.
  inputTokensPerSecond: 32_000,
  outputTokensPerSecond: 8_000,
  capacityMultiplier: 1, // Scales BOTH simultaneously attainable service rates.
  inputPeakShare: 0.5,
  outputPeakShare: 0.5,
  inputAnnualGrowth: 0,
  outputAnnualGrowth: 0,
  opexAnnualGrowth: 0, // Escalates fixed OpEx and extra local per-token cost rates.
  perMillionLocalInputCost: 0, // Only costs excluded from monthlyOpex.
  perMillionLocalOutputCost: 0,
  apiCacheHit: 0.75,
  fallbackCacheHit: null, // Null follows apiCacheHit; otherwise independently priced.
  peakHours: 154, // Synthetic 22 non-holiday workdays × 7 h, not a calendar.
  offHours: 566,
  availability: 0.98,
  headroom: 0.8, // Fraction retained for scheduling: 0.8 leaves 20% spare.
  capex: 1_600_000, // Illustrative 1.5M equipment + 0.1M implementation budget.
  monthlyOpex: 24_000, // Includes the example's energy, labor and maintenance.
  horizon: 36,
  annualDiscount: 0.1,
  terminalResidual: 240_000,
  apiAnnualChange: 0,
  prices: Object.freeze({
    peak: Object.freeze({ hit: 0.04, miss: 2, output: 8 }),
    off: Object.freeze({ hit: 0.02, miss: 1, output: 4 }),
  }),
});

export const examplePreset = Object.freeze({
  id: 'deepseek-v4.1-flash-example',
  metadata: Object.freeze({
    presetId: 'deepseek-v4.1-flash-example',
    name: 'DeepSeek V4.1 Flash', version: 'DeepSeek-V4.1-Flash-1017', currency: 'CNY',
    checkedAt: '2026-10-20',
    pricingSource: 'https://api-docs.deepseek.com/zh-cn/quick_start/pricing/',
    status: 'Illustrative budget and service-capacity assumptions; not a quote or hardware benchmark',
  }),
  parameters: exampleParameters,
});

// Kept for dated regression reproduction; the public calculator imports defaults.js.
export {customDefaults} from './defaults.js';
