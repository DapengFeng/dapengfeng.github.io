import sharp from 'sharp';
import QRCode from 'qrcode';
import path from 'node:path';
import {site,categories,escape as e} from './config.mjs';
import {sharePaths} from './sharing.mjs';
import {art} from './templates.mjs';

// Rasterize locally at build time. No QR/image service receives a visitor's URL.
async function textLayer(text,width,size,color='#eff3e9',bold=false){
 return sharp({text:{text:`<span foreground="${color}">${e(text)}</span>`,font:`Noto Sans CJK SC ${bold?'Bold ':''}${size}`,width,rgba:true,wrap:'word-char',spacing:8}}).png().toBuffer({resolveWithObject:true});
}
async function fitText(text,width,size,height,color,bold){
 let result;
 do{result=await textLayer(text,width,size,color,bold);size-=2;}while(result.info.height>height&&size>=20);
 if(result.info.height>height)throw Error(`Sharing text exceeds poster bounds (${result.info.height}px): ${text}`);
 return result;
}
function excerpt(text,length){const chars=Array.from(text);return chars.length<=length?text:chars.slice(0,length).join('').replace(/\s+\S*$/,'')+'…';}
export async function renderShareAssets(post,write){
 const paths=sharePaths(post),url=site.url+post.url;
 const qr=await QRCode.toBuffer(url,{type:'png',errorCorrectionLevel:'M',margin:4,scale:8,color:{dark:'#142019',light:'#ffffff'}});
 await write(paths.qr.slice(1),qr);
 const category=categories.find(c=>c.id===post.category);
 let visual;
 if(post.cover){
  if(!post.cover.startsWith('/assets/content/'))throw Error('Share cover must use a local content asset');
  const root=path.resolve('content/assets'),file=path.resolve(root,post.cover.slice('/assets/content/'.length));
  if(!file.startsWith(root+path.sep))throw Error('Invalid share cover path');
  visual=await sharp(file).resize(804,310,{fit:'cover'}).png().toBuffer();
 }else{
  visual=await sharp(Buffer.from(art(post.art,post.slug+'-share').replace('<svg ', '<svg xmlns="http://www.w3.org/2000/svg" width="804" height="310" '))).resize(804,310).png().toBuffer();
 }
 const qrSmall=await sharp(qr).resize(206,206,{kernel:'nearest'}).png().toBuffer();
 for(const mode of ['en','zh','both']){
  const title=mode==='zh'?post.title:mode==='en'?post.titleEn:post.titleEn+'\n'+post.title;
  const description=(scale)=>mode==='zh'?excerpt(post.description,Math.floor(180*scale)):mode==='en'?excerpt(post.descriptionEn,Math.floor(340*scale)):excerpt(post.descriptionEn,Math.floor(190*scale))+'\n'+excerpt(post.description,Math.floor(85*scale));
  const heading=await fitText(title,804,mode==='both'?42:50,228,'#f0f3eb',true);
  let detail,summaryScale=1;
  do{detail=await textLayer(description(summaryScale),804,28,'#c1cbbb');summaryScale*=.85;}while(detail.info.height>917-533-heading.info.height&&summaryScale>.15);
  if(detail.info.height>917-533-heading.info.height)throw Error(`Summary exceeds poster bounds: ${post.slug}/${mode}`);
  const name=mode==='en'?site.author:mode==='zh'?site.authorZh:site.author+' / '+site.authorZh;
  const subtitle=mode==='en'?category.en:mode==='zh'?category.name:category.en+' / '+category.name;
  const brand=await textLayer('FENG / '+(mode==='zh'?'知识实验室':mode==='en'?'KNOWLEDGE LAB':'KNOWLEDGE LAB / 知识实验室'),650,25,'#c0f47b');
  const categoryText=await textLayer(post.date+'   ·   '+subtitle,804,23,'#c1cbbb');
  const author=await textLayer(name,500,29);
  const scan=await textLayer(mode==='en'?'Scan to read':mode==='zh'?'扫码阅读全文':'Scan to read / 扫码阅读全文',500,25,'#c0f47b');
  const domain=await textLayer(new URL(site.url).hostname,500,23,'#c1cbbb');
  const base=Buffer.from('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="1200"><rect width="900" height="1200" fill="#101810"/><path d="M48 98H852M48 929H852" stroke="#47553f"/><circle cx="838" cy="54" r="7" fill="#c0f47b"/></svg>');
  const layers=[{input:brand.data,left:48,top:36},{input:visual,left:48,top:124},{input:categoryText.data,left:48,top:461},{input:heading.data,left:48,top:511},{input:detail.data,left:48,top:533+heading.info.height},{input:author.data,left:48,top:970},{input:scan.data,left:48,top:1024},{input:domain.data,left:48,top:1110},{input:qrSmall,left:646,top:958}];
  await write(paths.poster(mode).slice(1),await sharp(base).composite(layers).png().toBuffer());
 }
}
