import {watch} from 'node:fs';
import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
import {serve} from './serve.mjs';
const run=promisify(execFile);
async function build(){const result=await run(process.execPath,['scripts/build.mjs']);process.stdout.write(result.stdout);}
await build();serve();
let timer,busy=false,again=false;
async function rebuild(){if(busy){again=true;return;}busy=true;try{await build();console.log('Updated. Refresh your browser.');}catch(e){console.error(e.stderr||e);}finally{busy=false;if(again){again=false;rebuild();}}}
for(const dir of ['src','content','scripts'])watch(dir,{recursive:true},()=>{clearTimeout(timer);timer=setTimeout(rebuild,180);});
