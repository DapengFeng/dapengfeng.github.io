// One fictional full year, CNY 10,000. All opening balances are zero.
export function statements(collected = 80) {
  if (!Number.isFinite(collected) || collected < 0 || collected > 100) throw new RangeError('Collection must be between 0 and 100.');
  const revenue = 100, purchases = 50, supplierCash = 35, costOfSales = 30;
  const cashOpex = 20, depreciation = 6, interest = 4, tax = 10;
  const contributed = 100, debt = 40, capex = 60;
  const grossProfit = revenue - costOfSales;
  const operatingProfit = grossProfit - cashOpex - depreciation;
  const profit = operatingProfit - interest - tax;
  const receivables = revenue - collected, inventory = purchases - costOfSales;
  const payables = purchases - supplierCash, equipment = capex - depreciation;
  const operatingCash = collected - supplierCash - cashOpex - interest - tax;
  const investingCash = -capex, financingCash = contributed + debt;
  const cash = operatingCash + investingCash + financingCash;
  const equity = contributed + profit;
  return {collected, revenue, grossProfit, operatingProfit, profit, receivables, inventory,
    payables, equipment, debt, cash, equity, operatingCash, investingCash, financingCash,
    assets: cash + receivables + inventory + equipment,
    liabilitiesAndEquity: payables + debt + equity,
    indirectCash: profit + depreciation - receivables - inventory + payables};
}
