(() => {
  const $ = s => document.querySelector(s);
  const $$ = s => [...document.querySelectorAll(s)];
  const language = $('#site-language');
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const bi = (en, zh) => `<span class="i18n"><span data-lang="en">${esc(en)}</span><span data-lang="zh">${esc(zh)}</span></span>`;
  window.LabI18n = {bi, esc, get language(){return document.documentElement.dataset.language||'both';}};
  function applyLanguage(value, manual=false) {
    document.documentElement.dataset.language=value;
    document.documentElement.lang=value==='zh'?'zh-CN':'en';
    language.querySelectorAll('[data-language-choice]').forEach(button=>button.setAttribute('aria-pressed',String(button.dataset.languageChoice===value)));
    if(manual){document.documentElement.dataset.languageSource='manual';try { localStorage.setItem('feng-language',value); } catch {}}
    $$('[data-en][data-zh]').forEach(el=>el.textContent=value==='both'?`${el.dataset.en} / ${el.dataset.zh}`:el.dataset[value]);
    for(const attr of ['placeholder','aria-label','title']) $$(`[data-${attr}-en]`).forEach(el=>el.setAttribute(attr,value==='both'?`${el.getAttribute(`data-${attr}-en`)} / ${el.getAttribute(`data-${attr}-zh`)}`:el.getAttribute(`data-${attr}-${value}`)));
    window.dispatchEvent(new CustomEvent('languagechange',{detail:value}));
  }
  applyLanguage(document.documentElement.dataset.language||'both');
  let stopLocation;
  language.addEventListener('click',event=>{const button=event.target.closest('[data-language-choice]');if(button){stopLocation?.();applyLanguage(button.dataset.languageChoice,true);}});
  // A manual choice wins; automatic results last for this tab's session only.
  if(document.documentElement.dataset.languageSource==='browser'){
    const fallback=document.documentElement.dataset.language,controller=new AbortController();
    stopLocation=()=>controller.abort();
    const timeout=setTimeout(stopLocation,2500);
    fetch('https://api.country.is/',{signal:controller.signal,credentials:'omit',referrerPolicy:'no-referrer',cache:'no-store'})
      .then(response=>{if(!response.ok)throw Error('Country lookup unavailable');return response.json();})
      .then(data=>{if(typeof data.country!=='string'||! /^[A-Z]{2}$/.test(data.country))throw Error('Invalid country');return ['CN','HK','MO','TW'].includes(data.country)?'zh':'en';})
      .catch(()=>fallback)
      .then(value=>{
        if(document.documentElement.dataset.languageSource==='manual')return;
        // Do not overwrite a preference chosen in another tab during the request.
        try{const saved=localStorage.getItem('feng-language');if(['en','zh','both'].includes(saved)){applyLanguage(saved);document.documentElement.dataset.languageSource='manual';return;}}catch{}
        try{sessionStorage.setItem('feng-auto-language',value);}catch{}
        document.documentElement.dataset.languageSource='auto';applyLanguage(value);
      }).finally(()=>{clearTimeout(timeout);stopLocation=null;});
  }
  $('.mobile-menu')?.addEventListener('click',e=>{const open=$('.desktop-nav').classList.toggle('open');e.currentTarget.setAttribute('aria-expanded',String(open));});
  const dialog=$('#knowledge-search'), input=$('#global-search'), results=$('.search-results');
  let searchData=null,loading=null;
  async function loadSearch(){if(searchData)return searchData;if(!loading)loading=fetch('/search-index.json').then(r=>{if(!r.ok)throw Error('Search unavailable');return r.json();}).then(data=>searchData=data).catch(error=>{loading=null;throw error;});return loading;}
  async function search(){
    const query=input.value.trim().toLowerCase();
    try { const all=await loadSearch(); if(query!==input.value.trim().toLowerCase())return;
      const tokens=query.split(/\s+/).filter(Boolean);
      const found=all.filter(p=>tokens.every(q=>[p.title,p.titleEn,p.description,p.descriptionEn,...p.tags,p.searchText].join(' ').toLowerCase().includes(q))).slice(0,12);
      results.innerHTML=found.length?found.map(p=>`<a class="search-result" href="${esc(p.url)}"><small>${esc(p.date)} · ${bi(p.categoryEn,p.categoryZh)}</small><strong>${bi(p.titleEn||p.title,p.title)}</strong><span>${bi(p.descriptionEn||p.description,p.description)}</span></a>`).join(''):`<div class="empty-state">${bi('No matching notes. Try another keyword.','没有找到匹配的笔记，试试其他关键词。')}</div>`;
    }catch{results.innerHTML=bi('Search could not load. Please try again.','搜索暂时无法加载，请重试。');}
  }
  function openSearch(){if(!dialog.open)dialog.showModal();input.focus();search();}
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
