(() => {
 const dialog=document.getElementById('article-support-dialog');
 if(!dialog||typeof dialog.showModal!=='function')return;
 const link=dialog.querySelector('.support-paypal-link');
 const custom=dialog.querySelector('[data-support-custom]');
 let customEdited=false;
 const selection=()=>dialog.querySelector('[name="support-amount"]:checked');
 const updateAmount=()=>{
  const choice=selection()?.value,isCustom=!choice;
  custom.closest('.support-custom').classList.toggle('is-selected',isCustom);
  custom.setCustomValidity('');
  const amount=Number(isCustom?custom.value:choice);
  const valid=!isCustom||(custom.checkValidity()&&amount>0&&Number.isSafeInteger(Math.round(amount*100)));
  if(valid){
   link.href=`${link.dataset.supportProfile}/${amount}USD`;
   link.removeAttribute('aria-disabled');
   custom.removeAttribute('aria-invalid');
  }else{
   const language=document.documentElement.dataset.language;
   const en='Enter a positive amount with up to two decimal places.',zh='请输入大于 0、最多两位小数的金额。';
   custom.setCustomValidity(language==='zh'?zh:language==='en'?en:`${en} / ${zh}`);
   link.removeAttribute('href');
   link.setAttribute('aria-disabled','true');
   custom.setAttribute('aria-invalid','true');
  }
 };
 for(const radio of dialog.querySelectorAll('[name="support-amount"]'))radio.addEventListener('change',updateAmount);
 const selectCustom=()=>{
  for(const radio of dialog.querySelectorAll('[name="support-amount"]'))radio.checked=false;
  updateAmount();
 };
 // A text field with a decimal keyboard supports consistent selection across browsers.
 // Selecting the initial text alone never changes which amount is payable.
 const selectInitialText=()=>{if(!customEdited)custom.select();};
 custom.addEventListener('focus',selectInitialText);
 custom.addEventListener('click',()=>{selectCustom();selectInitialText();});
 custom.addEventListener('input',()=>{customEdited=true;selectCustom();});
 link.addEventListener('click',event=>{
  updateAmount();
  if(!link.hasAttribute('href')){event.preventDefault();custom.reportValidity();custom.focus();}
 });
 custom.addEventListener('keydown',event=>{if(event.key==='Enter'){event.preventDefault();selectCustom();link.click();}});
 updateAmount();
 let opener;
 for(const button of document.querySelectorAll('[data-support-open]')){
  button.hidden=false;
  button.addEventListener('click',()=>{
   opener=button;
   dialog.showModal();
   (selection()||custom).focus({preventScroll:true});
  });
 }
 dialog.addEventListener('keydown',event=>{
  if(event.key!=='Tab')return;
  const focusable=[...dialog.querySelectorAll('button,input,a[href],[tabindex="0"]')]
   .filter(node=>!node.disabled&&node.tabIndex>=0&&node.getClientRects().length&&(node.type!=='radio'||node.checked));
  const first=focusable[0],last=focusable.at(-1);
  if(event.shiftKey&&document.activeElement===first){event.preventDefault();last.focus();}
  else if(!event.shiftKey&&document.activeElement===last){event.preventDefault();first.focus();}
 });
 dialog.querySelector('[data-support-close]').addEventListener('click',()=>dialog.close());
 // Require both ends of a click to be outside, so a drag out of the dialog does not dismiss it.
 const outside=event=>{
  const r=dialog.getBoundingClientRect();
  return event.target===dialog&&(event.clientX<r.left||event.clientX>r.right||event.clientY<r.top||event.clientY>r.bottom);
 };
 let backdrop=false;
 dialog.addEventListener('pointerdown',event=>{backdrop=outside(event);});
 dialog.addEventListener('click',event=>{if(backdrop&&outside(event))dialog.close();backdrop=false;});
 dialog.addEventListener('close',()=>opener?.focus({preventScroll:true}));
})();
