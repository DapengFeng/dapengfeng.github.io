// Native bridge to giscus's widget protocol. GitHub auth/publishing runs in giscus;
// the parent stores only the encrypted giscus session, as its official client does.
// https://github.com/giscus/giscus/blob/main/client.ts
(() => {
 const root=document.querySelector('#article-discussions');
 if(!root)return;
 const ORIGIN='https://giscus.app',SESSION='giscus-session';
 const shell=root.querySelector('.discussion-shell'),launcher=root.querySelector('.discussion-launcher'),close=root.querySelector('.discussion-close');
 const threads=[...root.querySelectorAll('.discussion-thread')],buttons=[...root.querySelectorAll('[data-discussion-kind]')];
 const frames=new Map(),timers=new Map();
 let current=threads[0],floating=false,footerVisible=false;
 const language=()=>document.documentElement.dataset.language||'both';
 const text=(en,zh)=>language()==='zh'?zh:language()==='en'?en:`${en} / ${zh}`;
 const widgetLanguage=()=>language()==='zh'?'zh-CN':'en';
 const storage={get(){try{return JSON.parse(localStorage.getItem(SESSION)||'""');}catch{return '';}},set(value){try{value?localStorage.setItem(SESSION,JSON.stringify(value)):localStorage.removeItem(SESSION);}catch{}}};
 const callback=new URL(location.href),returnedSession=callback.searchParams.get('giscus');
 if(returnedSession){storage.set(returnedSession);callback.searchParams.delete('giscus');history.replaceState(history.state,'',callback.pathname+callback.search+callback.hash);}
 let session=returnedSession||storage.get();
 if(typeof session!=='string')session='';
 const state=new Map();
 function fitViewport(){
  const viewport=window.visualViewport;
  if(!viewport)return;
  shell.style.setProperty('--discussion-viewport-height',`${viewport.height}px`);
  shell.style.setProperty('--discussion-keyboard-inset',`${Math.max(0,innerHeight-viewport.height-viewport.offsetTop)}px`);
 }
 window.visualViewport?.addEventListener('resize',fitViewport);
 window.visualViewport?.addEventListener('scroll',fitViewport);
 fitViewport();
 function status(thread,kind){
  state.set(thread,kind);
  const node=thread.querySelector('.discussion-status');
  node.hidden=kind==='ready';
  const frame=frames.get(thread);
  if(frame)frame.hidden=kind==='setup';
  node.textContent=kind==='loading'?text('Loading GitHub discussion…','正在加载 GitHub 讨论…'):kind==='setup'?text('Comments are waiting for the site owner to connect GitHub.','评论功能等待站点作者完成 GitHub 接入。'):kind==='empty'?text('No posts yet. Start this article’s thread below.','暂无发言，可在下方发起这篇文章的讨论。'):text('Unable to load the discussion. Retry or view it on GitHub.','暂时无法加载讨论，请重试或前往 GitHub 查看。');
  thread.querySelector('.discussion-retry').hidden=!['error','setup'].includes(kind);
 }
 function stopTimer(thread){clearTimeout(timers.get(thread));timers.delete(thread);}
 function send(frame,config){frame.contentWindow?.postMessage({giscus:{setConfig:config}},ORIGIN);}
 function load(thread){
  if(frames.has(thread))return;
  status(thread,'loading');
  const origin=new URL(location.pathname,location.origin);
  origin.searchParams.set('discussion',thread.dataset.kind);
  origin.hash='article-discussions';
  const params=new URLSearchParams({origin:origin.href,session,repo:root.dataset.repo,repoId:root.dataset.repoId,category:thread.dataset.category,categoryId:thread.dataset.categoryId,term:thread.dataset.term,strict:'1',description:root.dataset.description,backLink:root.dataset.backlink,theme:'dark_dimmed',reactionsEnabled:'1',emitMetadata:'1',inputPosition:'top'});
  const frame=document.createElement('iframe');
  frame.title=text('GitHub '+thread.dataset.kind,'GitHub '+({comment:'评论',idea:'想法',discussion:'讨论'}[thread.dataset.kind]));
  frame.className='discussion-frame';frame.allow='clipboard-write';frame.referrerPolicy='strict-origin-when-cross-origin';
  frame.src=`${ORIGIN}/${widgetLanguage()}/widget?${params}`;
  frames.set(thread,frame);
  thread.querySelector('.discussion-embed').append(frame);
  // A load event can mean an error page. Only protocol messages establish readiness.
  timers.set(thread,setTimeout(()=>status(thread,'error'),20000));
  frame.addEventListener('load',()=>send(frame,{lang:widgetLanguage()}));
 }
 function select(thread){
  current=thread;
  threads.forEach(t=>t.hidden=t!==thread);
  buttons.forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.discussionKind===thread.dataset.kind)));
  load(thread);
 }
 function toggle(open){
  floating=open;
  if(open){root.style.minHeight=`${root.offsetHeight}px`;shell.classList.add('is-floating');shell.setAttribute('role','region');shell.setAttribute('aria-labelledby','discussions-title');}
  else {shell.classList.remove('is-floating');shell.removeAttribute('role');shell.removeAttribute('aria-labelledby');root.style.minHeight='';}
  close.hidden=!open;launcher.setAttribute('aria-expanded',String(open));
  launcher.hidden=open||footerVisible;
  if(open){load(current);close.focus({preventScroll:true});}
  else if(!launcher.hidden)launcher.focus({preventScroll:true});
  else buttons.find(b=>b.getAttribute('aria-pressed')==='true').focus({preventScroll:true});
 }
 buttons.forEach(b=>b.addEventListener('click',()=>select(threads.find(t=>t.dataset.kind===b.dataset.discussionKind))));
 launcher.addEventListener('click',()=>toggle(true));close.addEventListener('click',()=>toggle(false));
 document.addEventListener('keydown',event=>{if(event.key==='Escape'&&floating)toggle(false);});
 threads.forEach(thread=>thread.querySelector('.discussion-retry').addEventListener('click',()=>{
  // Reload is explicit and only offered after an error; changing config to the
  // same term does not invalidate giscus's query cache.
  const frame=frames.get(thread);
  if(frame){status(thread,'loading');stopTimer(thread);frame.src=frame.src;timers.set(thread,setTimeout(()=>status(thread,'error'),20000));}
  else load(thread);
 }));
 window.addEventListener('message',event=>{
  if(event.origin!==ORIGIN||!event.data||typeof event.data!=='object'||!event.data.giscus)return;
  const thread=threads.find(t=>frames.get(t)?.contentWindow===event.source);
  if(!thread)return;
  const data=event.data.giscus,frame=frames.get(thread);
  if(typeof data.resizeHeight==='number'&&Number.isFinite(data.resizeHeight)&&data.resizeHeight>0){
   frame.style.height=`${Math.min(100000,Math.max(220,data.resizeHeight))}px`;
   if(state.get(thread)==='loading'){stopTimer(thread);status(thread,'ready');}
  }
  if(data.discussion){
   stopTimer(thread);status(thread,'ready');
  }
  if(data.signOut===true){
   session='';storage.set('');
   for(const f of frames.values()){const url=new URL(f.src);url.searchParams.delete('session');f.src=url.href;}
  }
  if(typeof data.error==='string'){
   stopTimer(thread);
   if(/Bad credentials|Invalid state value|State has expired/.test(data.error)&&session){
    session='';storage.set('');
    for(const f of frames.values()){const url=new URL(f.src);url.searchParams.delete('session');f.src=url.href;}
   }
   status(thread,data.error.includes('Discussion not found')?'empty':/not installed|installation/i.test(data.error)?'setup':'error');
  }
 });
 window.addEventListener('languagechange',()=>{
  for(const thread of threads){if(state.has(thread))status(thread,state.get(thread));const f=frames.get(thread);if(f){send(f,{lang:widgetLanguage()});f.title=text('GitHub '+thread.dataset.kind,'GitHub '+({comment:'评论',idea:'想法',discussion:'讨论'}[thread.dataset.kind]));}}
 });
 const observer=new IntersectionObserver(entries=>{
  footerVisible=entries[0].isIntersecting;
  launcher.hidden=floating||footerVisible;
  if(footerVisible)load(current);
 },{rootMargin:'0px'});
 observer.observe(root);
 const selected=threads.find(t=>t.dataset.kind===callback.searchParams.get('discussion'));
 if(returnedSession||location.hash==='#article-discussions'){if(selected)current=selected;select(current);toggle(true);}
 else launcher.hidden=false;
})();
