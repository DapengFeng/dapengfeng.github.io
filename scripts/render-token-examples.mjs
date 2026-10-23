import assert from 'node:assert/strict';
import {readFile, writeFile} from 'node:fs/promises';
import {teachingExample} from '../content/assets/llm-economics/token-examples.js';
import {tokenExamplesHTML} from '../content/assets/llm-economics/token-examples-view.js';

const article = new URL('../content/posts/llm-deployment-roi.html', import.meta.url);
const source = await readFile(article, 'utf8');
const start = '<!-- token-examples:start -->';
const end = '<!-- token-examples:end -->';
assert.equal(source.split(start).length, 2, 'exactly one teaching figure start');
assert.equal(source.split(end).length, 2, 'exactly one teaching figure end');
const begin = source.indexOf(start) + start.length;
const finish = source.indexOf(end);
assert.ok(finish > begin);
const generated = `\n${tokenExamplesHTML(teachingExample, {mode: 'teaching'})}\n`;
if (process.argv.includes('--check')) {
  assert.equal(source.slice(begin, finish), generated, 'regenerate the static token comparison');
  console.log('Static token comparison matches the calculation and live view.');
} else {
  await writeFile(article, source.slice(0, begin) + generated + source.slice(finish));
  console.log('Updated the static token comparison in the article source.');
}
