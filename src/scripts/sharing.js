(() => {
 const dialog=document.getElementById('article-share-dialog');
 if(!dialog||typeof dialog.showModal!=='function')return;
 const triggers=[...document.querySelectorAll('[data-share-open]')];
 const field=dialog.querySelector('#article-share-text'),status=dialog.querySelector('[data-share-status]');
 const download=dialog.querySelector('[data-share-download]'),native=dialog.querySelector('[data-share-native]'),fileButton=dialog.querySelector('[data-share-file]');
 const manual=dialog.querySelector('.share-manual'),manualField=dialog.querySelector('#share-manual-text');
 const url=dialog.dataset.shareUrl,drafts=new Map();
 let mode,opener,platform='',preparedFile,controller,requestId=0,busy=false;
 const lang=()=>['en','zh'].includes(document.documentElement.dataset.language)?document.documentElement.dataset.language:'both';
 const text=(en,zh)=>mode==='en'?en:mode==='zh'?zh:`${en} / ${zh}`;
 const title=()=>text(dialog.dataset.titleEn,dialog.dataset.titleZh);
 const summary=()=>text(dialog.dataset.summaryEn,dialog.dataset.summaryZh);
 function message(en='',zh=''){status.textContent=text(en,zh);}
 function updateImage(image,path){
  image.alt=text(image.dataset.altEn,image.dataset.altZh);
  if(path)image.dataset.src=path;
  // Nothing loads until the reader opens the corresponding sharing view.
  if(dialog.open&&!image.closest('[hidden]'))image.src=image.dataset.src;
 }
 function update(){
  mode=lang();preparedFile=undefined;fileButton.hidden=true;
  field.value=drafts.get(mode)??`${title()}\n\n${summary()}\n\n${url}`;
  download.href=`${dialog.dataset.posterBase}-${mode}.png`;
  download.download=`${dialog.dataset.posterBase.split('/').pop()}-${mode}.png`;
  // X counts CJK characters more heavily; leave room for the URL.
  let short='',weight=0;
  for(const char of title()){const size=char.codePointAt(0)>0x10ff?2:1;if(weight+size>215){short+='…';break;}short+=char;weight+=size;}
  dialog.querySelector('[data-share-platform=x]').href='https://x.com/intent/tweet?'+new URLSearchParams({url,text:short});
  dialog.querySelector('[data-share-platform=telegram]').href='https://t.me/share/url?'+new URLSearchParams({url,text:title()});
  updateImage(dialog.querySelector('[data-share-poster] img'),download.getAttribute('href'));
  updateImage(dialog.querySelector('[data-share-qr] img'));
  manual.hidden=true;status.textContent='';
  if(dialog.open)prepareFile();
 }
 async function prepareFile(){
  const id=++requestId;controller?.abort();preparedFile=undefined;fileButton.hidden=true;
  if(!navigator.share||!navigator.canShare)return;
  controller=new AbortController();
  try{
   const response=await fetch(download.href,{signal:controller.signal});
   if(!response.ok)throw Error('Image unavailable');
   const blob=await response.blob();if(blob.type!=='image/png')throw Error('Invalid image');
   const file=new File([blob],download.download,{type:'image/png'});
   if(id!==requestId||!dialog.open)return;
   if(navigator.canShare({files:[file]})){preparedFile=file;fileButton.hidden=false;}
  }catch{/* Link, text, QR, and direct download remain available. */}
 }
 function showPlatform(next){
  platform=next;
  field.hidden=platform==='wechat';
  dialog.querySelector('.share-text-label').hidden=platform==='wechat';
  dialog.querySelectorAll('button[data-share-platform]').forEach(button=>button.setAttribute('aria-pressed',String(button.dataset.sharePlatform===platform)));
  dialog.querySelector('[data-share-media]').hidden=!platform;
  dialog.querySelector('[data-share-qr]').hidden=platform!=='wechat';
  dialog.querySelector('[data-share-poster]').hidden=platform!=='rednote';
  if(platform)updateImage(dialog.querySelector(platform==='wechat'?'[data-share-qr] img':'[data-share-poster] img'));
 }
 triggers.forEach(button=>{
  button.hidden=false;
  button.addEventListener('click',()=>{
   opener=button;update();showPlatform('');dialog.showModal();prepareFile();
  });
 });
 dialog.querySelector('[data-share-close]').addEventListener('click',()=>dialog.close());
 dialog.addEventListener('click',event=>{if(event.target===dialog){const box=dialog.getBoundingClientRect();if(event.clientX<box.left||event.clientX>box.right||event.clientY<box.top||event.clientY>box.bottom)dialog.close();}});
 dialog.addEventListener('close',()=>{controller?.abort();requestId++;preparedFile=undefined;fileButton.hidden=true;opener?.focus({preventScroll:true});});
 dialog.querySelectorAll('button[data-share-platform]').forEach(button=>button.addEventListener('click',()=>showPlatform(button.dataset.sharePlatform)));
 field.addEventListener('input',()=>drafts.set(mode,field.value));
 dialog.querySelectorAll('[data-share-action]').forEach(button=>button.addEventListener('click',async()=>{
  const value=button.dataset.shareAction==='copy-link'?url:field.value;
  button.disabled=true;
  try{
   if(!navigator.clipboard?.writeText)throw Error('Clipboard unavailable');
   await navigator.clipboard.writeText(value);manual.hidden=true;message('Copied','已复制');
  }catch{
   manual.hidden=false;manualField.value=value;manualField.focus();manualField.select();
   message('Select the text below to copy.','请选中下方文字复制。');
  }finally{button.disabled=false;}
 }));
 native.hidden=typeof navigator.share!=='function';
 function share(withFile){
  if(busy)return;
  if(withFile&&!preparedFile)return;
  busy=true;native.disabled=true;fileButton.disabled=true;
  const data=withFile?{files:[preparedFile],title:title(),text:field.value}:{title:title(),text:summary(),url};
  // Call synchronously within the click: no fetch or await before navigator.share.
  let task;
  try{task=navigator.share(data);}catch(error){task=Promise.reject(error);}
  Promise.resolve(task).then(()=>message('Handed to system sharing','已交给系统分享')).catch(error=>{
   if(error.name==='AbortError'){status.textContent='';return;}
   message('Sharing unavailable. Copy the link or save the image.','暂时无法调用系统分享，请复制链接或保存图片。');
  }).finally(()=>{busy=false;native.disabled=false;fileButton.disabled=false;});
 }
 native.addEventListener('click',()=>share(false));fileButton.addEventListener('click',()=>share(true));
 dialog.querySelectorAll('.share-media img').forEach(image=>image.addEventListener('error',()=>message('Image unavailable. You can still copy the link.','图片暂时无法加载，仍可复制链接。')));
 window.addEventListener('languagechange',update);update();
})();
