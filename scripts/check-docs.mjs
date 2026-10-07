import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

// Check repository documentation without installing the site's/browser's dependencies.
export async function checkDocs(root=process.cwd()){
 const files=['README.md'];
 try{
  await fs.access(path.join(root,'AGENTS.md'));
  files.push('AGENTS.md');
 }catch(error){if(error.code!=='ENOENT')throw error;}
 async function walk(directory){
  for(const item of await fs.readdir(path.join(root,directory),{withFileTypes:true})){
   const file=path.posix.join(directory,item.name);
   if(item.isDirectory())await walk(file);
   else if(file.endsWith('.md'))files.push(file);
  }
 }
 await walk('docs');
 const errors=[];
 for(const file of files){
  const text=new TextDecoder('utf-8',{fatal:true}).decode(await fs.readFile(path.join(root,file)));
  if(!text.trim())errors.push(`${file}: empty document`);
  let fence=null;
  for(const [i,line]of text.split('\n').entries()){
   const where=`${file}:${i+1}`;
   if(/^(<{7}|={7}|>{7})(\s|$)/.test(line))errors.push(`${where}: unresolved merge conflict`);
   const marker=/^\s{0,3}(`{3,}|~{3,})/.exec(line)?.[1];
   if(marker){
    if(!fence)fence=marker;
    else if(marker[0]===fence[0]&&marker.length>=fence.length)fence=null;
    continue;
   }
   if(fence)continue;
   // Inline links and reference definitions; remote URLs and fragment-only links are excluded.
   const links=[...line.matchAll(/\]\(\s*(?:<([^>]+)>|([^\s)]+))/g)].map(m=>m[1]||m[2]);
   const reference=/^\s{0,3}\[[^\]]+\]:\s*(?:<([^>]+)>|(\S+))/.exec(line);
   if(reference)links.push(reference[1]||reference[2]);
   for(const link of links){
    if(/^(?:[a-z][a-z\d+.-]*:|\/\/|#)/i.test(link))continue;
    try{
     const pathname=decodeURIComponent(link.split(/[?#]/)[0]);
     const target=pathname.startsWith('/')?path.join(root,pathname):path.resolve(root,path.dirname(file),pathname);
     await fs.access(target);
    }catch{errors.push(`${where}: missing local link target ${link}`);}
   }
  }
  if(fence)errors.push(`${file}: unclosed code fence`);
 }
 if(errors.length)throw Error(errors.join('\n'));
 return files.length;
}

if(process.argv[1]===fileURLToPath(import.meta.url))console.log(`Documentation checks passed (${await checkDocs()} files).`);
