(() => {
  const $ = s => document.querySelector(s);
  const $$ = s => [...document.querySelectorAll(s)];
  const monthlyVolume = $('[data-current-month]');
  if (monthlyVolume) {
    // Match the site's calendar even when visitors browse from another time zone.
    const month = new Intl.DateTimeFormat('en', {month:'numeric', timeZone:'Asia/Shanghai'});
    const updateVolume = () => {
      const value = `VOL. ${month.format(new Date()).padStart(3,'0')}`;
      if (monthlyVolume.textContent !== value) monthlyVolume.textContent = value;
    };
    updateVolume();
    setInterval(updateVolume, 60000);
    window.addEventListener('pageshow', updateVolume);
    document.addEventListener('visibilitychange', () => { if (!document.hidden) updateVolume(); });
  }
  const language = $('#site-language');
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const bi = (en, zh) => `<span class="i18n"><span data-lang="en" lang="en">${esc(en)}</span><span data-lang="zh" lang="zh-CN">${esc(zh)}</span></span>`;
  window.LabI18n = {bi, esc, get language(){return document.documentElement.dataset.language||'both';}};
  function applyLanguage(value, manual=false) {
    document.documentElement.dataset.language=value;
    document.documentElement.lang=value==='zh'?'zh-CN':'en';
    language?.querySelectorAll('[data-language-choice]').forEach(button=>button.setAttribute('aria-pressed',String(button.dataset.languageChoice===value)));
    if(manual){document.documentElement.dataset.languageSource='manual';try { localStorage.setItem('feng-language',value); } catch {}}
    $$('[data-en][data-zh]').forEach(el=>el.textContent=value==='both'?`${el.dataset.en} / ${el.dataset.zh}`:el.dataset[value]);
    for(const attr of ['placeholder','aria-label','title']) $$(`[data-${attr}-en]`).forEach(el=>el.setAttribute(attr,value==='both'?`${el.getAttribute(`data-${attr}-en`)} / ${el.getAttribute(`data-${attr}-zh`)}`:el.getAttribute(`data-${attr}-${value}`)));
    $$('[data-reading-minutes]').forEach(el=>{
      const minutes=el.dataset[value==='both'?'minutesBoth':value==='zh'?'minutesZh':'minutesEn'];
      if(!minutes)return;
      el.innerHTML=bi(`${minutes} min read`,`${minutes} 分钟阅读`);
    });
    window.dispatchEvent(new CustomEvent('languagechange',{detail:value}));
  }
  applyLanguage(document.documentElement.dataset.language||'both');
  language?.addEventListener('click',event=>{const button=event.target.closest('[data-language-choice]');if(button)applyLanguage(button.dataset.languageChoice,true);});
  $('.mobile-menu')?.addEventListener('click',e=>{const open=$('.desktop-nav').classList.toggle('open');e.currentTarget.setAttribute('aria-expanded',String(open));});
  const dialog=$('#knowledge-search'), input=$('#global-search'), results=$('.search-results');
  const header=$('.lab-header'),navigation=$('.desktop-nav'),menu=$('.mobile-menu');
  let previousScroll=Math.max(0,scrollY),scrollDirection=0,scrollDistance=0,headerFrame=0;
  function revealHeader(){
    header?.classList.remove('is-retracted');
    previousScroll=Math.max(0,scrollY);scrollDistance=0;scrollDirection=0;
  }
  if(header){
    // Translate the sticky header without changing its height or shifting the page.
    const updateHeader=()=>{
      headerFrame=0;
      const y=Math.max(0,Math.min(scrollY,document.documentElement.scrollHeight-innerHeight));
      const delta=y-previousScroll;previousScroll=y;
      const keyboardFocus=header.contains(document.activeElement)&&document.activeElement.matches(':focus-visible');
      if(y<=header.offsetHeight||navigation?.classList.contains('open')||dialog?.open||keyboardFocus){revealHeader();return;}
      if(!delta)return;
      const direction=Math.sign(delta);
      scrollDistance=direction===scrollDirection?scrollDistance+Math.abs(delta):Math.abs(delta);
      scrollDirection=direction;
      if(direction>0&&scrollDistance>=36)header.classList.add('is-retracted');
      if(direction<0&&scrollDistance>=12)header.classList.remove('is-retracted');
    };
    window.addEventListener('scroll',()=>{if(!headerFrame)headerFrame=requestAnimationFrame(updateHeader);},{passive:true});
    header.addEventListener('focusin',revealHeader);
    menu?.addEventListener('click',revealHeader);
    navigation?.addEventListener('click',event=>{
      if(event.target.closest('a')){navigation.classList.remove('open');menu?.setAttribute('aria-expanded','false');revealHeader();}
    });
    document.addEventListener('keydown',event=>{
      if(event.key==='Tab')revealHeader();
      if(event.key==='Escape'&&navigation?.classList.contains('open')){
        navigation.classList.remove('open');menu?.setAttribute('aria-expanded','false');menu?.focus();revealHeader();
      }
    });
    dialog?.addEventListener('close',revealHeader);
    window.addEventListener('pageshow',revealHeader);
    new ResizeObserver(()=>{
      document.documentElement.style.setProperty('--site-header-height',`${header.offsetHeight}px`);
      revealHeader();
    }).observe(header);
  }
  let searchData=null,loading=null;
  async function loadSearch(){if(searchData)return searchData;if(!loading)loading=fetch('/search-index.json').then(r=>{if(!r.ok)throw Error('Search unavailable');return r.json();}).then(data=>searchData=data).catch(error=>{loading=null;throw error;});return loading;}
  function searchScore(record,tokens){
    const fields=[
      [[record.title,record.titleEn].join(' '),100],
      [(record.tags||[]).join(' '),60],
      [[record.description,record.descriptionEn].join(' '),20],
      [[record.searchText,record.searchTextEn,record.searchTextZh].join(' '),1],
    ].map(([text,weight])=>[text.toLowerCase(),weight]);
    let score=0;
    for(const token of tokens){
      const match=fields.find(([text])=>text.includes(token));
      if(!match)return -1;
      score+=match[1];
    }
    return score;
  }
  function excerpt(record,locale,tokens){
    const fallback=locale==='en'?(record.descriptionEn||record.description):(record.description||record.descriptionEn);
    // Avoid presenting a concatenation of both editions as a search excerpt.
    const body=record[locale==='en'?'searchTextEn':'searchTextZh'];
    if(!body||!tokens.length)return fallback||'';
    const normalized=body.toLowerCase();
    const matches=tokens.map(token=>normalized.indexOf(token)).filter(index=>index>=0);
    if(!matches.length)return fallback||'';
    const position=Math.min(...matches),length=locale==='zh'?100:180;
    let start=Math.max(0,position-Math.floor(length/4)),end=Math.min(body.length,start+length);
    if(locale==='en'){
      const boundary=body.lastIndexOf(' ',start);if(boundary>=0&&start-boundary<24)start=boundary+1;
      const tail=body.indexOf(' ',end);if(tail>=0&&tail-end<24)end=tail;
    }
    return `${start?'…':''}${body.slice(start,end).trim()}${end<body.length?'…':''}`;
  }
  async function search(){
    const query=input.value.trim().toLowerCase();
    try { const all=await loadSearch(); if(query!==input.value.trim().toLowerCase())return;
      const tokens=query.split(/\s+/).filter(Boolean);
      const found=all.map((record,order)=>({record,order,score:searchScore(record,tokens)})).filter(item=>item.score>=0).sort((a,b)=>b.score-a.score||a.order-b.order).slice(0,12).map(item=>item.record);
      results.innerHTML=found.length?found.map(p=>`<a class="search-result" href="${esc(p.url)}"><small>${p.date?esc(p.date)+' · ':''}${bi(p.categoryEn||'Note',p.categoryZh||'笔记')}</small><strong>${bi(p.titleEn||p.title,p.title||p.titleEn)}</strong><span>${bi(excerpt(p,'en',tokens),excerpt(p,'zh',tokens))}</span></a>`).join(''):`<div class="empty-state">${bi('No matching notes. Try another keyword.','没有找到匹配的笔记，试试其他关键词。')}</div>`;
    }catch{results.innerHTML=bi('Search could not load. Please try again.','搜索暂时无法加载，请重试。');}
  }
  function openSearch(){revealHeader();if(!dialog.open)dialog.showModal();input.focus();search();}
  $$('.search-trigger').forEach(b=>b.addEventListener('click',openSearch));
  $('[data-close-search]')?.addEventListener('click',()=>dialog.close());
  dialog?.addEventListener('click',event=>{if(event.target===dialog){const r=dialog.getBoundingClientRect();if(event.clientX<r.left||event.clientX>r.right||event.clientY<r.top||event.clientY>r.bottom)dialog.close();}});
  let searchTimer;input?.addEventListener('input',()=>{clearTimeout(searchTimer);searchTimer=setTimeout(search,100);});
  document.addEventListener('keydown',e=>{if(e.key==='Escape'&&dialog.open){e.preventDefault();dialog.close();return;}if((e.metaKey||e.ctrlKey)&&e.key.toLowerCase()==='k'){e.preventDefault();dialog.open?dialog.close():openSearch();}});
  const grid=$('[data-library]');
  if(grid){
    const cards=[...grid.children], categoryButtons=$$('[data-category-filter]'), queryInput=$('#library-search'), sort=$('#article-sort');
    const params=new URLSearchParams(location.search);let category=params.get('category')||'all';if(!categoryButtons.some(b=>b.dataset.categoryFilter===category))category='all';
    if(queryInput)queryInput.value=params.get('q')||'';
    function filter(sync=false){
      const query=(queryInput?.value||'').trim().toLowerCase(),tokens=query.split(/\s+/).filter(Boolean);
      let count=0;for(const card of cards){const shown=(category==='all'||card.dataset.category===category)&&tokens.every(t=>card.dataset.search.includes(t));card.hidden=!shown;if(shown)count++;}
      categoryButtons.forEach(b=>{const active=b.dataset.categoryFilter===category;b.classList.toggle('active',active);b.setAttribute('aria-pressed',String(active));});
      const chinese=document.documentElement.dataset.language==='zh';
      const titleKey=chinese?'titleZh':'titleEn',collator=new Intl.Collator(chinese?'zh-CN':'en',{numeric:true,sensitivity:'base'});
      cards.sort((a,b)=>sort.value==='title'?collator.compare(a.dataset[titleKey],b.dataset[titleKey]):sort.value==='oldest'?a.dataset.date.localeCompare(b.dataset.date):b.dataset.date.localeCompare(a.dataset.date)).forEach(c=>grid.append(c));
      $('.empty-state').hidden=count>0;$('#result-count').innerHTML=bi(`${count} notes in the collection`,`共 ${count} 篇笔记`);
      if(sync){const url=new URL(location.href);category==='all'?url.searchParams.delete('category'):url.searchParams.set('category',category);query?url.searchParams.set('q',query):url.searchParams.delete('q');history.replaceState(null,'',url);}
    }
    categoryButtons.forEach(b=>b.addEventListener('click',()=>{category=b.dataset.categoryFilter;filter(true);}));queryInput?.addEventListener('input',()=>filter(true));sort?.addEventListener('change',()=>filter());
    window.addEventListener('languagechange',()=>filter());
    $('[data-reset-filters]')?.addEventListener('click',()=>{category='all';if(queryInput)queryInput.value='';filter(true);});filter();
  }
  // Retire the previous Jekyll cache worker, if this browser used the old site.
  if('serviceWorker' in navigator)navigator.serviceWorker.getRegistrations().then(regs=>regs.filter(r=>new URL(r.active?.scriptURL||location.href).pathname==='/sw.js'||new URL(r.active?.scriptURL||location.href).pathname==='/assets/scripts/sw.js').forEach(r=>r.unregister())).catch(()=>{});
})();
// A quiet pointer light follows the cards, without changing their layout.
if(matchMedia('(pointer: fine)').matches&&!matchMedia('(prefers-reduced-motion: reduce)').matches){
 document.querySelectorAll('.knowledge-card,.category-tile').forEach(card=>card.addEventListener('pointermove',event=>{const r=card.getBoundingClientRect();card.style.setProperty('--pointer-x',`${event.clientX-r.left}px`);card.style.setProperty('--pointer-y',`${event.clientY-r.top}px`);}));
}
