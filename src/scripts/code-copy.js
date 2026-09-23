(() => {
 const content=document.getElementById('article-content');if(!content)return;
 const mounted=new WeakSet();
 const pair=(en,zh)=>`<span data-lang="en" lang="en">${en}</span><span data-lang="zh" lang="zh-CN">${zh}</span>`;
 const icon='<svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><g class="copy-symbol"><rect x="8" y="8" width="12" height="12" rx="2"></rect><path d="M16 8V5a2 2 0 0 0-2-2H5a2 2 0 0 0-2 2v9a2 2 0 0 0 2 2h3"></path></g><path class="copied-symbol" d="m5 12 4 4L19 6"></path></svg>';
 function textOf(pre){
  const copy=pre.cloneNode(true);
  copy.querySelectorAll('.line-number,[aria-hidden="true"]').forEach(n=>n.remove());
  const language=document.documentElement.dataset.language;
  if(language==='en'||language==='zh')copy.querySelectorAll(`[data-lang]:not([data-lang="${language}"])`).forEach(n=>n.remove());
  const rows=[...copy.querySelectorAll('.code-line')];
  return rows.length?rows.map(n=>n.textContent).join('\n'):copy.textContent;
 }
 function mount(target,source){
  if(mounted.has(target))return;mounted.add(target);
  const block=document.createElement('div');block.className='code-copy-block';
  const tools=document.createElement('div');tools.className='code-copy-tools';
  const status=document.createElement('span');status.className='code-copy-status sr-only';status.setAttribute('role','status');
  const button=document.createElement('button');button.type='button';button.className='code-copy';
  button.innerHTML=icon+'<span class="code-copy-label sr-only"></span>';
  tools.append(status,button);target.before(block);block.append(tools,target);
  let timer;
  function feedback(copied=false){
   button.dataset.copied=String(copied);
   const en=copied?'Copied':'Copy code',zh=copied?'已复制':'复制代码';
   button.querySelector('.code-copy-label').innerHTML=pair(en,zh);
   const language=document.documentElement.dataset.language;
   button.title=language==='en'?en:language==='zh'?zh:`${en} / ${zh}`;
  }
  feedback();window.addEventListener('languagechange',()=>feedback(button.dataset.copied==='true'));
  button.addEventListener('click',async()=>{
   clearTimeout(timer);button.disabled=true;const text=source();
   try{
    if(!navigator.clipboard?.writeText)throw Error('Clipboard unavailable');
    await navigator.clipboard.writeText(text);
    block.querySelector('.code-copy-manual')?.remove();feedback(true);
    status.innerHTML=pair('Code copied.','代码已复制。');
    timer=setTimeout(()=>{feedback();status.textContent='';},2000);
   }catch{
    feedback();status.innerHTML=pair('Automatic copying is unavailable. Select the code below to copy it manually.','无法自动复制，请选中下方代码手动复制。');
    let manual=block.querySelector('.code-copy-manual');
    if(!manual){
     manual=document.createElement('label');manual.className='code-copy-manual';
     manual.innerHTML=pair('Select and copy code','选中并复制代码');
     const field=document.createElement('textarea');field.readOnly=true;field.rows=5;field.spellcheck=false;
     manual.append(field);block.append(manual);
    }
    const field=manual.querySelector('textarea');field.value=text;field.focus();field.select();
   }finally{button.disabled=false;}
  });
 }
 function scan(){
  // Old article-specific controls are superseded by the shared icon.
  content.querySelectorAll('button[data-copy]').forEach(n=>n.hidden=true);
  content.querySelectorAll('pre:not(.compiler-highlight)').forEach(pre=>{
   if(pre.hidden||pre.closest('[aria-hidden="true"]'))return;
   mount(pre,()=>textOf(pre));
   if(!pre.hasAttribute('tabindex'))pre.tabIndex=0;
  });
  content.querySelectorAll('.compiler-editor').forEach(editor=>mount(editor.closest('.compiler-editor-surface')||editor,()=>editor.value));
 }
 scan();
 // Also cover results and code blocks inserted by an article's interactive examples.
 new MutationObserver(records=>{
  if(records.some(record=>[...record.addedNodes].some(n=>n.nodeType===1&&(n.matches('pre,.compiler-editor')||n.querySelector('pre,.compiler-editor')))))scan();
 }).observe(content,{childList:true,subtree:true});
})();
