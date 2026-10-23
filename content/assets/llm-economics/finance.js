/**
 * Currency-neutral investment arithmetic. All amounts must use the same currency.
 * Month zero pays initialInvestment. Each supplied net flow occurs at month end;
 * terminalValue is an additional net receipt (or disposal cost) in the last month.
 * Revenue, replacement equipment and other costs belong in monthlyNetCashFlows
 * when applicable. This module makes no assumption about models, tokens or vendors.
 */
function paybackHistory(months, field) {
  const crossing = months.findIndex(row => row[field] >= 0);
  const reversalMonths = [];
  for (let i = Math.max(1, crossing + 1); crossing >= 0 && i < months.length; i++) {
    if (months[i - 1][field] >= 0 && months[i][field] < 0) reversalMonths.push(months[i].month);
  }
  let sustainedMonth = null;
  if (months.at(-1)[field] >= 0) {
    let i = months.length - 1;
    while (i > 0 && months[i - 1][field] >= 0) i--;
    sustainedMonth = months[i].month;
  }
  return { month: crossing < 0 ? null : months[crossing].month, reversalMonths, sustainedMonth };
}

/**
 * Payback excludes terminalValue and reports later reversals. capitalROI means
 * undiscounted net gain / initialInvestment, not annualized return or return on
 * all later capital outlays. With initialInvestment=0 it is undefined (null).
 * Investment limits hold the supplied flows and terminal value fixed.
 */
export function evaluateCashFlows({ initialInvestment, monthlyNetCashFlows, terminalValue = 0, annualDiscount = 0 } = {}) {
  if (!Number.isFinite(initialInvestment) || initialInvestment < 0) throw new RangeError('initialInvestment must be finite and nonnegative');
  if (!Array.isArray(monthlyNetCashFlows) || monthlyNetCashFlows.length < 1 || !monthlyNetCashFlows.every(Number.isFinite)) {
    throw new RangeError('monthlyNetCashFlows must be a nonempty array of finite net amounts');
  }
  if (!Number.isFinite(terminalValue)) throw new RangeError('terminalValue must be finite');
  if (!Number.isFinite(annualDiscount) || annualDiscount <= -1) throw new RangeError('annualDiscount must be finite and greater than -1');

  const horizon = monthlyNetCashFlows.length;
  const months = [{
    month: 0, netCashFlow: -initialInvestment, terminalValue: 0, discountFactor: 1,
    cashFlow: -initialInvestment, discountedCashFlow: -initialInvestment,
    cumulativePreTerminalNet: -initialInvestment, cumulativeDiscountedPreTerminalNet: -initialInvestment,
    cumulativeNet: -initialInvestment, cumulativeDiscountedNet: -initialInvestment,
  }];
  let presentValueMonthlyNetCashFlows = 0;
  for (let month = 1; month <= horizon; month++) {
    const netCashFlow = monthlyNetCashFlows[month - 1];
    const discountFactor = (1 + annualDiscount) ** (-month / 12);
    const residual = month === horizon ? terminalValue : 0;
    const previous = months.at(-1);
    const cumulativePreTerminalNet = previous.cumulativePreTerminalNet + netCashFlow;
    const cumulativeDiscountedPreTerminalNet = previous.cumulativeDiscountedPreTerminalNet + netCashFlow * discountFactor;
    const row = {
      month, netCashFlow, terminalValue: residual, discountFactor,
      cashFlow: netCashFlow + residual, discountedCashFlow: (netCashFlow + residual) * discountFactor,
      cumulativePreTerminalNet, cumulativeDiscountedPreTerminalNet,
      cumulativeNet: cumulativePreTerminalNet + residual,
      cumulativeDiscountedNet: cumulativeDiscountedPreTerminalNet + residual * discountFactor,
    };
    if (!Object.values(row).every(Number.isFinite)) throw new RangeError('Projected cash flows or discounted values are not finite');
    months.push(row);
    presentValueMonthlyNetCashFlows += netCashFlow * discountFactor;
  }
  const presentValueTerminalValue = terminalValue * (1 + annualDiscount) ** (-horizon / 12);
  const simple = paybackHistory(months, 'cumulativePreTerminalNet');
  const discounted = paybackHistory(months, 'cumulativeDiscountedPreTerminalNet');
  const final = months.at(-1);
  return {
    horizon, initialInvestment, annualDiscount,
    monthlyDiscount: (1 + annualDiscount) ** (1 / 12) - 1, months,
    presentValueMonthlyNetCashFlows, presentValueTerminalValue,
    maxInitialInvestment: presentValueMonthlyNetCashFlows + presentValueTerminalValue,
    maxInitialInvestmentWithoutTerminal: presentValueMonthlyNetCashFlows,
    preTerminalNetCash: final.cumulativePreTerminalNet,
    terminalNetCash: final.cumulativeNet,
    NPV: final.cumulativeDiscountedNet, preTerminalNPV: final.cumulativeDiscountedPreTerminalNet,
    capitalROI: initialInvestment > 0 ? final.cumulativeNet / initialInvestment : null,
    simplePaybackMonth: simple.month, discountedPaybackMonth: discounted.month,
    simplePaybackReversalMonths: simple.reversalMonths, discountedPaybackReversalMonths: discounted.reversalMonths,
    sustainedSimplePaybackMonth: simple.sustainedMonth, sustainedDiscountedPaybackMonth: discounted.sustainedMonth,
  };
}
