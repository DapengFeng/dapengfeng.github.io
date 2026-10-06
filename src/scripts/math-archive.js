import {dateInfo} from './daily-math.js';
const root=document.querySelector('[data-math-archive]');
if(root){
 let timer;
 function update(){
  clearTimeout(timer);const {date,nextChangeAt}=dateInfo();let count=0;
  root.querySelectorAll('[data-entry-date]').forEach(row=>{row.hidden=row.dataset.entryDate>date||(!root.hasAttribute('data-year-archive')&&count>=30);if(!row.hidden)count++;});
  root.querySelectorAll('[data-math-year]').forEach(link=>{link.hidden=link.dataset.mathYear>date.slice(0,4);});
  root.querySelector('[data-math-archive-empty]').hidden=count>0;
  timer=setTimeout(update,Math.min(3600000,Math.max(100,nextChangeAt-Date.now()+30)));
 }
 window.addEventListener('pageshow',update);document.addEventListener('visibilitychange',()=>{if(!document.hidden)update();});window.addEventListener('pagehide',()=>clearTimeout(timer));update();
}
