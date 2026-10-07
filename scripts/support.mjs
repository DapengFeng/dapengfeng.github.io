import {escape as e} from './config.mjs';
import {pair} from './i18n.mjs';

export function validateSupport(input){
 const url=input?.paypal?.url??'';
 // A public PayPal.Me profile; do not accept arbitrary redirects, query strings or checkout amounts.
 if(typeof url!=='string'||(url&&(url!==url.trim()||!/^https:\/\/paypal\.me\/[a-zA-Z0-9_-]+\/?$/.test(url))))throw Error('PayPal support URL must be an HTTPS paypal.me profile link');
 return {paypal:{url}};
}
export function loadSupport(env=process.env){return validateSupport({paypal:{url:env.PAYPAL_ME_URL?.trim()??''}});}

const svg=body=>`<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${body}</svg>`;
const icons={
 coffee:svg('<path d="M4 8h13v8a4 4 0 0 1-4 4H8a4 4 0 0 1-4-4Zm13 1h1a3 3 0 0 1 0 6h-1M3 23h16M7 2v2m4-3v3m4-2v2"/>'),
 close:svg('<path d="m6 6 12 12M6 18 18 6"/>'),
 paypal:'<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M9 5h7c4 0 5 2 4 6-1 4-3 5-7 5h-2l-1 6H5Z" fill="#74c7ff"/><path d="M6 2h7c4 0 6 2 5 6-1 4-3 5-7 5H9l-1 6H3Z" fill="#367fca"/><path d="m9 5-2 8h4c4 0 6-1 7-5l.3-2C17.7 5.3 16.6 5 15 5Z" fill="#205b9f"/></svg>'
};
const label=(en,zh)=>`aria-label="${e(en+' / '+zh)}" data-aria-label-en="${e(en)}" data-aria-label-zh="${e(zh)}" title="${e(en+' / '+zh)}" data-title-en="${e(en)}" data-title-zh="${e(zh)}"`;
export function supportTrigger(){
 return `<button type="button" class="article-support-trigger" data-support-open aria-haspopup="dialog" aria-controls="article-support-dialog" ${label('Buy me a coffee','请我喝杯咖啡')} hidden>${icons.coffee}</button>`;
}
export function supportDialog(config){
 if(!config.paypal.url)return '';
 const profile=config.paypal.url.replace(/\/$/,'');
 return `<dialog id="article-support-dialog" class="article-support-dialog" aria-labelledby="article-support-title"><div class="support-heading"><h2 id="article-support-title">${pair('Buy me a coffee','请我喝杯咖啡')}</h2><button type="button" data-support-close ${label('Close support panel','关闭赞赏面板')}>${icons.close}</button></div><div class="support-paypal">${icons.paypal}<span class="support-recipient">${e(profile.split('/')[3])}</span>
<fieldset class="support-amounts"><legend class="sr-only">${pair('Amount in US dollars','金额（美元）')}</legend>${[1,3,5].map(amount=>`<label><input class="sr-only" type="radio" name="support-amount" value="${amount}"${amount===3?' checked':''}><span>US$${amount}</span></label>`).join('')}<label class="support-custom"><span aria-hidden="true">US$</span><input type="text" data-support-custom value="10" pattern="(?:[0-9]+(?:[.][0-9]{1,2})?|[.][0-9]{1,2})" inputmode="decimal" autocomplete="off" spellcheck="false" required ${label('Custom amount in US dollars','自定义金额（美元）')}></label></fieldset>
<a class="support-paypal-link" data-support-profile="${e(profile)}" href="${e(profile)}/3USD" target="_blank" rel="noopener noreferrer" tabindex="0" ${label('Continue to PayPal','前往 PayPal')}>${pair('Continue to PayPal','前往 PayPal',true)}</a></div></dialog><noscript><p class="support-fallback"><a href="${e(config.paypal.url)}" target="_blank" rel="noopener noreferrer">${pair('Support via PayPal','通过 PayPal 支持')}</a></p></noscript>`;
}
