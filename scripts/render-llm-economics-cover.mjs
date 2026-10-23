import fs from 'node:fs/promises';
import sharp from 'sharp';
import {evaluateCashFlows} from '../content/assets/llm-economics/finance.js';

// A generic cash-flow illustration, independent of any model or API quote.
const x = month => 76 + month / 36 * 1048;
const y = cash => 484 - (cash + 100000) / 180000 * 220;
const curves = [1500, 2800, 4500].map((saving, index) => {
  const r = evaluateCashFlows({
    initialInvestment: 100000, monthlyNetCashFlows: Array(36).fill(saving),
    annualDiscount: .1,
  });
  return `<path d="${r.months.map((row, i) => `${i ? 'L' : 'M'}${x(i).toFixed(2)},${y(row.cumulativeDiscountedPreTerminalNet).toFixed(2)}`).join(' ')}" fill="none" stroke="${['#efbd88', '#94d4ec', '#bcebc5'][index]}" stroke-width="4"/>`;
}).join('');
const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="630" viewBox="0 0 1200 630">
<rect width="1200" height="630" fill="#101b24"/>
<g font-family="DejaVu Sans, sans-serif">
<text x="72" y="62" fill="#94d4ec" font-size="18" letter-spacing="3">MODEL DEPLOYMENT · A COST COMPARISON</text>
<text x="72" y="139" fill="#f3f1e9" font-size="49" font-weight="600">API OR YOUR OWN HARDWARE?</text>
<text x="72" y="207" fill="#f3f1e9" font-size="36">MONTHLY SAVINGS · PAYBACK · NPV</text>
<path d="M76 ${y(0)}H1124" stroke="#667e8b" stroke-dasharray="5 8"/>
<path d="M76 264V501H1124" fill="none" stroke="#34464f"/>
${curves}
<text x="76" y="${y(0) - 14}" fill="#b9c8d1" font-size="17">BREAK-EVEN</text>
<text x="76" y="529" fill="#b9c8d1" font-size="18">0</text>
<text x="1124" y="529" fill="#b9c8d1" font-size="18" text-anchor="end">36 MONTHS</text>
<text x="72" y="586" fill="#b9c8d1" font-size="20">Model prices · Your workload · Measured capacity</text>
</g></svg>`;
const out = new URL('../content/assets/llm-economics/cover.png', import.meta.url);
await fs.mkdir(new URL('.', out), {recursive: true});
await sharp(Buffer.from(svg)).png().toFile(out.pathname);
console.log('Rendered generic economics cover (illustrative cash flows).');
