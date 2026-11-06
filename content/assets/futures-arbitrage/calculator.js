/**
 * Close an equal-quantity, equal-multiplier calendar spread at one shared time.
 * Prices use the same currency and unit; cost is the complete round-trip cost
 * for both legs in that currency. This models P&L, not margin or daily cash needs.
 */
export function calculateCalendarSpread({
  nearOpen,
  farOpen,
  nearClose,
  farClose,
  multiplier = 10,
  lots = 1,
  cost = 40,
  longNear = true,
} = {}) {
  for (const [name, value] of Object.entries({
    nearOpen, farOpen, nearClose, farClose, multiplier, lots, cost,
  })) {
    if (typeof value !== 'number' || !Number.isFinite(value)) {
      throw new RangeError(`${name} must be a finite number`);
    }
  }
  if (multiplier <= 0) throw new RangeError('multiplier must be positive');
  if (!Number.isInteger(lots) || lots <= 0) {
    throw new RangeError('lots must be a positive integer');
  }
  if (cost < 0) throw new RangeError('cost must be nonnegative');
  if (typeof longNear !== 'boolean') throw new RangeError('longNear must be a boolean');

  const quantity = multiplier * lots;
  const direction = longNear ? 1 : -1;
  const spreadOpen = farOpen - nearOpen;
  const spreadClose = farClose - nearClose;
  const nearPnl = direction * (nearClose - nearOpen) * quantity;
  const farPnl = -direction * (farClose - farOpen) * quantity;
  const gross = nearPnl + farPnl;
  const net = gross - cost;
  const breakEvenFar = nearClose + spreadOpen - direction * cost / quantity;

  return { nearPnl, farPnl, gross, net, spreadOpen, spreadClose, breakEvenFar, quantity };
}
