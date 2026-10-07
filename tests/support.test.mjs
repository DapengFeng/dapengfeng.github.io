import test from 'node:test';
import assert from 'node:assert/strict';
import {load} from 'cheerio';
import {loadSupport,validateSupport,supportDialog} from '../scripts/support.mjs';

test('PayPal accepts only a public HTTPS profile, not arbitrary payment redirects',()=>{
 for(const url of ['http://paypal.me/person','https://paypal.me.evil.invalid/person','https://paypal.me@evil.invalid/person','javascript:alert(1)','https://paypal.me/person/100','https://paypal.me/person?redirect=elsewhere','https://paypal.me/person#fragment','https://paypal.me/person\n']){
  assert.throws(()=>validateSupport({paypal:{url}}),/PayPal/);
 }
 assert.equal(validateSupport({paypal:{url:'https://paypal.me/example'}}).paypal.url,'https://paypal.me/example');
});

test('support reads the optional build variable and never falls back to a recipient',()=>{
 for(const env of [{},{PAYPAL_ME_URL:''},{PAYPAL_ME_URL:'  '}]){
  assert.equal(loadSupport(env).paypal.url,'');
  assert.equal(supportDialog(loadSupport(env)),'');
 }
 assert.equal(loadSupport({PAYPAL_ME_URL:' https://paypal.me/example\n'}).paypal.url,'https://paypal.me/example');
 assert.throws(()=>loadSupport({PAYPAL_ME_URL:'https://example.com/checkout'}),/PayPal/);
});

test('production exposes only PayPal, with no payment-method tabs or placeholders',()=>{
 const config=loadSupport({PAYPAL_ME_URL:'https://paypal.me/example'}),$=load(supportDialog(config));
 assert.equal($('[role=tablist],.support-unavailable,[data-support-image]').length,0);
 assert.doesNotMatch($('#article-support-dialog').text(),/WeChat|Alipay|微信|支付宝|暂未开通/);
 assert.equal($('.support-paypal-link').attr('href'),config.paypal.url+'/3USD');
 assert.equal($('.support-paypal-link').attr('rel'),'noopener noreferrer');
 assert.equal(supportDialog(validateSupport({paypal:{url:''}})),'');
 assert.match($('noscript').html(),/paypal\.me\/example/);
});
