import {test} from 'node:test';
import assert from 'node:assert/strict';
import vm from 'node:vm';
import fs from 'node:fs/promises';
import {load} from 'cheerio';
const context = {window:{}, document:{querySelectorAll:()=>[]}};
vm.runInNewContext(await fs.readFile('src/scripts/syntax.js','utf8'), context);
const {highlight} = context.window.CodeSyntax;
test('syntax coloring preserves exact source, including markup and unfinished edits', () => {
  for (const language of ['c++','rust','python','unknown']) {
    for (const source of ['<img src=x onerror=alert(1)> & "quoted"\n\t中文 😀', '/* unfinished\ncomment', '"unfinished', "'a\n", '', 'return 1;\n']) {
      const $ = load(`<code>${highlight(source,language)}</code>`);
      assert.equal($('code').text(),source);
      assert.equal($('img, script').length,0);
    }
  }
});
test('comments, raw strings and Python triple quotes keep their lexical boundaries', () => {
  for (const [language,source,kind] of [
    ['c++','/* int fake() {} */','comment'],
    ['rust','/* outer /* nested */ still a comment */','comment'],
    ['c++','R"tag(<img> "quoted"\n)tag"','string'],
    ['rust','r##"<script> "quoted"\n"##','string'],
    ['python','f"""text\nreturn False\n"""','string'],
    ['python','# return False','comment']
  ]) {
    const $=load(`<code>${highlight(source,language)}</code>`);
    assert.equal($('code').text(),source);
    assert.equal($('span').length,1);
    assert.equal($('span').attr('class'),`syntax-${kind}`);
  }
});
