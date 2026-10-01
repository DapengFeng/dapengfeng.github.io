// Fetch the local renderer only near this article's figure. Other pages pay no cost.
for(const root of document.querySelectorAll('[data-eye-viewer]')){
 const poster=root.querySelector('.eye-poster');
 const localizePoster=()=>{const lang=document.documentElement.dataset.language||'both';poster.alt=lang==='both'?`${poster.dataset.altEn} / ${poster.dataset.altZh}`:poster.getAttribute(`data-alt-${lang}`);};
 localizePoster();window.addEventListener('languagechange',localizePoster);
 let requested=false;
 const observer=new IntersectionObserver(async entries=>{
  if(requested||!entries.some(e=>e.isIntersecting))return;requested=true;observer.disconnect();
  try{
   const {mountEye}=await import('./eye-renderer.js');
   const viewer=await mountEye(root);
   if(viewer)root.eyeViewer=viewer;
   else root.dataset.eyeState='fallback';
  }catch{root.dataset.eyeState='fallback';}
 },{rootMargin:'350px'});
 observer.observe(root);
}
