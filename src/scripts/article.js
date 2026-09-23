(() => {
 const content=document.getElementById('article-content'),links=[...document.querySelectorAll('.article-toc nav a')],bar=document.querySelector('.article-progress');if(!content)return;
 let scheduled=false;
 function update(){
  scheduled=false;const top=content.getBoundingClientRect().top+scrollY,max=Math.max(1,content.offsetHeight-innerHeight),ratio=Math.min(1,Math.max(0,(scrollY-top)/max));
  bar.style.width=`${ratio*100}%`;document.getElementById('read-percentage').textContent=`READING · ${Math.round(ratio*100)}%`;
  const visible=links.filter(a=>a.getClientRects().length);let active=visible[0];for(const link of visible){const node=document.getElementById(decodeURIComponent(link.hash.slice(1)));if(node?.getBoundingClientRect().top<160)active=link;}
  links.forEach(a=>{a.classList.toggle('active',a===active);if(a===active)a.setAttribute('aria-current','location');else a.removeAttribute('aria-current');});
 }
 const queue=()=>{if(!scheduled){scheduled=true;requestAnimationFrame(update);}};
 window.addEventListener('scroll',queue,{passive:true});window.addEventListener('resize',queue);window.addEventListener('languagechange',()=>{document.querySelectorAll('.article-edition').forEach(e=>{if(e.getClientRects().length)e.style.marginTop='';});queue();});new ResizeObserver(queue).observe(content);
 // Formula source is retained at build time; never copy rendered MathML or glyphs.
 const mathPair=(en,zh)=>`<span class="i18n"><span data-lang="en" lang="en">${en}</span><span data-lang="zh" lang="zh-CN">${zh}</span></span>`;
 content.querySelectorAll('.formula-block').forEach(block=>{
  const button=block.querySelector('.formula-copy'),status=block.querySelector('.formula-status');
  let timer;
  const label=button.querySelector('.formula-copy-label');
  function feedback(copied=false){
   button.dataset.copied=String(copied);
   const en=copied?'Copied':'Copy LaTeX',zh=copied?'已复制':'复制 LaTeX';
   label.innerHTML=mathPair(en,zh);
   const language=document.documentElement.dataset.language;
   button.title=language==='en'?en:language==='zh'?zh:`${en} / ${zh}`;
  }
  feedback();
  window.addEventListener('languagechange',()=>feedback(button.dataset.copied==='true'));
  button.addEventListener('click',async()=>{
   clearTimeout(timer);button.disabled=true;
   try{
    if(!navigator.clipboard?.writeText)throw Error('Clipboard unavailable');
    await navigator.clipboard.writeText(block.dataset.latex);
    block.querySelector('.formula-manual')?.remove();
    feedback(true);status.innerHTML=mathPair('LaTeX copied.','已复制 LaTeX。');
    timer=setTimeout(()=>{feedback();status.textContent='';},2000);
   }catch{
    feedback();
    status.innerHTML=mathPair('Automatic copying is unavailable. Select the LaTeX below to copy it manually.','无法自动复制，请选中下方 LaTeX 手动复制。');
    let manual=block.querySelector('.formula-manual');
    if(!manual){
     manual=document.createElement('div');manual.className='formula-manual';
     const label=document.createElement('label');label.innerHTML=mathPair('Select and copy LaTeX','选中并复制 LaTeX');
     const source=document.createElement('textarea');source.readOnly=true;source.rows=3;source.setAttribute('aria-label','LaTeX');source.spellcheck=false;
     source.value=block.dataset.latex;label.append(source);manual.append(label);block.append(manual);
    }
    const source=manual.querySelector('textarea');source.focus();source.select();
   }finally{button.disabled=false;}
  });
 });
 document.getElementById('print-article').addEventListener('click',()=>print());
 function reveal(){if(!location.hash)return;let target;try{target=document.getElementById(decodeURIComponent(location.hash.slice(1)));}catch{return;}if(!target)return;for(let node=target.parentElement;node;node=node.parentElement)if(node.tagName==='DETAILS')node.open=true;}
 window.addEventListener('hashchange',reveal);reveal();update();
})();
