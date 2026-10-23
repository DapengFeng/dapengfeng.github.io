/** Compatibility entry point for the published DeepSeek example and its figures.
 * General callers use llm-economics/token-model.js with a complete explicit profile.
 */
import { evaluateTokenScenario, billTokens as billExplicitTokens } from '../llm-economics/token-model.js';
import { examplePreset } from '../llm-economics/presets.js';

export const defaults = examplePreset.parameters;
function exampleParameters(overrides = {}) {
  return {
    ...defaults, ...overrides,
    prices: {
      peak: { ...defaults.prices.peak, ...overrides.prices?.peak },
      off: { ...defaults.prices.off, ...overrides.prices?.off },
    },
  };
}
export function evaluate(overrides = {}) {
  return evaluateTokenScenario(exampleParameters(overrides));
}
export function billTokens(inputTokenCount, outputTokenCount, overrides = {}, bucket = 'peak', month = 1, route = 'baseline') {
  return billExplicitTokens(inputTokenCount, outputTokenCount, exampleParameters(overrides), bucket, month, route);
}
export default defaults;
