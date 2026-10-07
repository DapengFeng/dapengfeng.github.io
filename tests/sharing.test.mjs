import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import sharp from 'sharp';
import jsQR from 'jsqr';
import * as cheerio from 'cheerio';
import {articleSharing,sharePaths} from '../scripts/sharing.mjs';
import {site} from '../scripts/config.mjs';
const posts=JSON.parse(await fs.readFile('dist/search-index.json','utf8')).filter(p=>p.kind==='article');
async function decoded(file){const {data,info}=await sharp(file).ensureAlpha().raw().toBuffer({resolveWithObject:true});return jsQR(new Uint8ClampedArray(data),info.width,info.height)?.data;}
test('every article has local sharing assets and a QR code decoding to its canonical URL',async()=>{
 for(const p of posts){
  const $=cheerio.load(await fs.readFile('dist'+p.url,'utf8')),paths=sharePaths(p);
  assert.equal($('[data-share-open]').length,1,p.slug);
  assert.equal($('.article-action-dock [data-share-open]').length,1,p.slug);
  assert.equal($('#article-share-dialog').attr('data-share-url'),site.url+p.url);
  assert.equal(await decoded('dist'+paths.qr),site.url+p.url,p.slug);
  for(const mode of ['en','zh','both']){
   const info=await sharp('dist'+paths.poster(mode)).metadata();
   assert.equal(info.width,900);assert.equal(info.height,1200);assert.equal(info.format,'png');
  }
  assert.equal($('.share-media img[src]').length,0,'sharing images load only on demand');
 }
});
test('photo and technical posters retain scannable QR codes in every language',async()=>{
 for(const slug of ['jiuzhaigou-water-and-mountains','pytorch-01-what-is-pytorch']){
  const post=posts.find(p=>p.slug===slug);
  for(const mode of ['en','zh','both'])assert.equal(await decoded('dist'+sharePaths(post).poster(mode)),site.url+post.url);
 }
});
test('sharing metadata escapes markup and platform URLs preserve the canonical address',()=>{
 const post={slug:'a',url:'/blog/a.html',titleEn:'A "quote" <script>alert(1)</script>',title:'标题 & 内容',descriptionEn:'<img src=x onerror=alert(1)>',description:'摘要'};
 const $=cheerio.load(articleSharing(post));
 assert.equal($('script').length,0);assert.equal($('[onerror]').length,0);
 assert.equal($('#article-share-dialog').attr('data-title-en'),post.titleEn);
 const href=new URL($('[data-share-platform=x]').attr('href'));
 assert.equal(href.searchParams.get('url'),site.url+post.url);
 assert.equal(href.searchParams.get('text'),post.titleEn);
 $('a[target=_blank]').each((_,el)=>assert.match($(el).attr('rel'),/noopener/));
});
