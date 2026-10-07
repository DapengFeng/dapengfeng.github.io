// Native bridge to giscus's widget protocol. GitHub auth/publishing runs in giscus;
// the parent stores only the encrypted giscus session, as its official client does.
// https://github.com/giscus/giscus/blob/main/client.ts
(() => {
 const root=document.querySelector('#article-discussions');
 if(!root)return;
 const ORIGIN='https://giscus.app',SESSION='giscus-session';
 const shell=root.querySelector('.discussion-shell'),launcher=root.querySelector('.discussion-launcher'),close=root.querySelector('.discussion-close');
 const thread=root.querySelector('.discussion-thread');
 let frame,timer,frameLanguage,state,frameLoaded=false,floating=false,footerVisible=false;
 const language=()=>document.documentElement.dataset.language||'both';
 const text=(en,zh)=>language()==='zh'?zh:language()==='en'?en:`${en} / ${zh}`;
 const widgetLanguage=()=>language()==='zh'?'zh-CN':'en';
 const storage={get(){try{return JSON.parse(localStorage.getItem(SESSION)||'""');}catch{return '';}},set(value){try{value?localStorage.setItem(SESSION,JSON.stringify(value)):localStorage.removeItem(SESSION);}catch{}}};
 const callback=new URL(location.href),returnedSession=callback.searchParams.get('giscus');
 const legacyLink=['comment','idea','discussion'].includes(callback.searchParams.get('discussion'));
 if(returnedSession)storage.set(returnedSession);
 if(returnedSession||legacyLink){callback.searchParams.delete('giscus');callback.searchParams.delete('discussion');history.replaceState(history.state,'',callback.pathname+callback.search+callback.hash);}
 let session=returnedSession||storage.get();
 if(typeof session!=='string')session='';
 function fitViewport(){
  const viewport=window.visualViewport;
  if(!viewport)return;
  shell.style.setProperty('--discussion-viewport-height',`${viewport.height}px`);
  shell.style.setProperty('--discussion-keyboard-inset',`${Math.max(0,innerHeight-viewport.height-viewport.offsetTop)}px`);
 }
 window.visualViewport?.addEventListener('resize',fitViewport);
 window.visualViewport?.addEventListener('scroll',fitViewport);
 fitViewport();
 function status(kind){
  state=kind;
  const node=thread.querySelector('.discussion-status');
  node.hidden=kind==='ready'||kind==='empty';
  node.dataset.state=kind;
  if(frame)frame.hidden=kind==='setup';
  const message=kind==='loading'?text('Loading GitHub discussion…','正在加载 GitHub 讨论…'):kind==='setup'?text('Comments are waiting for the site owner to connect GitHub.','评论功能等待站点作者完成 GitHub 接入。'):kind==='empty'?text('No posts yet.','暂无发言。'):text('Unable to load the discussion. Please retry.','暂时无法加载讨论，请重试。');
  node.replaceChildren();
  const label=document.createElement('span');label.textContent=message;
  if(kind==='loading')label.className='sr-only';
  node.append(label);
  thread.querySelector('.discussion-retry').hidden=!['error','setup'].includes(kind);
 }
 function startLoading(){clearTimeout(timer);status('loading');timer=setTimeout(()=>status('error'),20000);}
 function syncLanguage(){
  // The initial URL already selects a locale; sending it again navigates twice.
  const lang=widgetLanguage();
  if(!frameLoaded||frameLanguage===lang)return;
  frameLanguage=lang;frame.contentWindow?.postMessage({giscus:{setConfig:{lang}}},ORIGIN);
 }
 function load(){
  if(frame)return;
  const origin=new URL(location.pathname,location.origin);
  origin.hash='article-discussions';
  const params=new URLSearchParams({origin:origin.href,session,repo:root.dataset.repo,repoId:root.dataset.repoId,category:thread.dataset.category,categoryId:thread.dataset.categoryId,term:thread.dataset.term,strict:'1',description:root.dataset.description,backLink:root.dataset.backlink,theme:new URL('/assets/giscus-theme.css',location.origin).href,reactionsEnabled:'1',emitMetadata:'1',inputPosition:'top'});
  frame=document.createElement('iframe');
  frame.title=text('GitHub discussion','GitHub 讨论');
  frame.className='discussion-frame';frame.allow='clipboard-write';frame.referrerPolicy='strict-origin-when-cross-origin';
  frameLanguage=widgetLanguage();
  frame.src=`${ORIGIN}/${frameLanguage}/widget?${params}`;
  frame.addEventListener('load',()=>{frameLoaded=true;syncLanguage();});
  startLoading();
  thread.querySelector('.discussion-embed').append(frame);
 }
 function reload(signOut=false){
  const url=new URL(frame.src);
  frameLanguage=widgetLanguage();url.pathname=`/${frameLanguage}/widget`;
  if(signOut)url.searchParams.delete('session');
  frameLoaded=false;startLoading();frame.src=url.href;
 }
 function toggle(open){
  floating=open;
  if(open){root.style.minHeight=`${root.offsetHeight}px`;shell.classList.add('is-floating');shell.setAttribute('role','region');shell.setAttribute('aria-labelledby','discussions-title');}
  else {shell.classList.remove('is-floating');shell.removeAttribute('role');shell.removeAttribute('aria-labelledby');root.style.minHeight='';}
  close.hidden=!open;launcher.setAttribute('aria-expanded',String(open));
  launcher.hidden=open||footerVisible;
  if(open){load();close.focus({preventScroll:true});}
  else if(!launcher.hidden)launcher.focus({preventScroll:true});
  else frame?.focus({preventScroll:true});
 }
 launcher.addEventListener('click',()=>toggle(true));close.addEventListener('click',()=>toggle(false));
 document.addEventListener('keydown',event=>{if(event.key==='Escape'&&floating)toggle(false);});
 thread.querySelector('.discussion-retry').addEventListener('click',()=>frame?reload():load());
 window.addEventListener('message',event=>{
  if(event.origin!==ORIGIN||!frame||event.source!==frame.contentWindow||!event.data||typeof event.data!=='object'||!event.data.giscus)return;
  const data=event.data.giscus;
  if(typeof data.resizeHeight==='number'&&Number.isFinite(data.resizeHeight)&&data.resizeHeight>0){
   frame.style.height=`${Math.min(100000,Math.max(220,data.resizeHeight))}px`;
   if(state==='loading'){clearTimeout(timer);status('ready');}
  }
  if(data.discussion){clearTimeout(timer);status('ready');}
  if(data.signOut===true){session='';storage.set('');reload(true);return;}
  if(typeof data.error==='string'){
   clearTimeout(timer);
   if(/Bad credentials|Invalid state value|State has expired/.test(data.error)&&session){
    session='';storage.set('');reload(true);return;
   }
   status(data.error.includes('Discussion not found')?'empty':/not installed|installation/i.test(data.error)?'setup':'error');
  }
 });
 window.addEventListener('languagechange',()=>{
  if(state)status(state);
  if(frame){syncLanguage();frame.title=text('GitHub discussion','GitHub 讨论');}
 });
 const observer=new IntersectionObserver(entries=>{
  footerVisible=entries[0].isIntersecting;
  launcher.hidden=floating||footerVisible;
  if(footerVisible)load();
 },{rootMargin:'0px'});
 observer.observe(root);
 // Warm the article's one editor without competing with the initial render.
 const speculativeAllowed=()=>!document.hidden&&!navigator.connection?.saveData&&!['slow-2g','2g'].includes(navigator.connection?.effectiveType);
 const warm=()=>{if(speculativeAllowed())load();};
 const nearby=new IntersectionObserver(entries=>{if(entries[0].isIntersecting)warm();},{rootMargin:'1500px 0px'});
 nearby.observe(root);
 launcher.addEventListener('pointerenter',warm);
 launcher.addEventListener('focus',warm);
 if(returnedSession||legacyLink||location.hash==='#article-discussions')toggle(true);
 else launcher.hidden=false;
})();
