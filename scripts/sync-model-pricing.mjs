#!/usr/bin/env node
import { readFile, writeFile, rename, unlink } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { updateCatalog, validateCatalog } from './lib/model-pricing.mjs';

const target = fileURLToPath(new URL('../content/assets/llm-economics/model-prices.json', import.meta.url));
const args = process.argv.slice(2);
if (args.some(arg => arg !== '--check')) throw new Error('Usage: node scripts/sync-model-pricing.mjs [--check]');
let previous;
try {
  previous = JSON.parse(await readFile(target, 'utf8'));
} catch (error) {
  if (error.code !== 'ENOENT') throw error;
}
if (args.includes('--check')) {
  validateCatalog(previous);
  console.log(`CC Switch upstream prices validated offline: ${previous.quotes.length} quotes from ${previous.coverage.importedProviders} providers.`);
} else {
  const { catalog, failures } = await updateCatalog(previous, async source => {
    const response = await fetch(source.url, {
      signal: AbortSignal.timeout(20_000),
      redirect: 'error',
      headers: { 'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0 (compatible; FENG-Knowledge-Lab-Price-Sync/2.0)' },
    });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const type = response.headers.get('content-type') ?? '';
    if (!/^application\/json\b/i.test(type)) throw new Error('Unexpected response content type');
    const declaredSize = Number(response.headers.get('content-length') ?? 0);
    if (declaredSize > 16_000_000) throw new Error('Response exceeds 16 MB limit');
    const chunks = [];
    let total = 0;
    for await (const chunk of response.body) {
      total += chunk.length;
      if (total > 16_000_000) throw new Error('Response exceeds 16 MB limit');
      chunks.push(chunk);
    }
    return JSON.parse(Buffer.concat(chunks).toString('utf8'));
  });
  const temporary = `${target}.${process.pid}.tmp`;
  try {
    await writeFile(temporary, `${JSON.stringify(catalog)}\n`, { flag: 'wx' });
    await rename(temporary, target);
  } finally {
    await unlink(temporary).catch(error => { if (error.code !== 'ENOENT') throw error; });
  }
  console.log(`CC Switch upstream prices updated: ${catalog.quotes.length} quotes from ${catalog.coverage.importedProviders} providers; ${catalog.sources.filter(s => s.status === 'ok').length}/${catalog.sources.length} sources refreshed.`);
  for (const failure of failures) console.warn(`${failure.sourceId}: retained its last verified snapshot (${failure.error}).`);
}
