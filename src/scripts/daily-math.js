// Shanghai calendar dates, independent of browser locale and daylight-saving time.
export const DAY_MS=86400000;
const OFFSET=8*3600000;
export function parseDate(date){
 if(typeof date!=='string'||!/^\d{4}-\d{2}-\d{2}$/.test(date))throw Error('Invalid date: '+date);
 const value=Date.parse(date+'T00:00:00Z');
 if(!Number.isFinite(value)||new Date(value).toISOString().slice(0,10)!==date)throw Error('Invalid date: '+date);
 return value;
}
export function addDays(date,days){return new Date(parseDate(date)+days*DAY_MS).toISOString().slice(0,10);}
export function dateInfo(timestamp=Date.now()){
 if(!Number.isFinite(timestamp))throw Error('Invalid date');
 const day=Math.floor((timestamp+OFFSET)/DAY_MS);
 return {date:new Date(day*DAY_MS).toISOString().slice(0,10),nextChangeAt:(day+1)*DAY_MS-OFFSET};
}
// Exact dates never wrap, including after the last scheduled entry or in a leap year.
export function dailySceneAt(timestamp,entries){
 const calendar=dateInfo(timestamp),index=entries.findIndex(entry=>entry.date===calendar.date);
 return {...calendar,index,scene:index<0?null:entries[index]};
}
export function fourier(t,terms=7){let value=0;for(let k=0;k<terms;k++){const n=2*k+1;value+=4/Math.PI*Math.sin(n*t)/n;}return value;}
export function lissajous(t){return [Math.sin(3*t),Math.cos(2*t)];}
export function standing(x,y,t){return Math.sin(2*x)*Math.sin(3*y)*Math.cos(Math.sqrt(13)*t);}
export function heatKernel(x,y,t){return Math.exp(-(x*x+y*y)/(4*t))/(4*Math.PI*t);}
export function mobius(u,v){return [(2+v*Math.cos(u/2))*Math.cos(u),(2+v*Math.cos(u/2))*Math.sin(u),v*Math.sin(u/2)];}
export function phyllotaxis(n){const angle=n*Math.PI*(3-Math.sqrt(5)),r=Math.sqrt(n);return [r*Math.cos(angle),r*Math.sin(angle)];}
