import http from 'node:http';
import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
export function serve(port=Number(process.env.PORT)||4173){
 const root=path.resolve('dist'),types={'.html':'text/html; charset=utf-8','.css':'text/css','.js':'text/javascript','.json':'application/json','.svg':'image/svg+xml','.png':'image/png','.jpg':'image/jpeg','.woff2':'font/woff2','.woff':'font/woff','.ttf':'font/ttf','.xml':'application/xml'};
 const server=http.createServer(async(req,res)=>{
  try{let pathname=decodeURIComponent(new URL(req.url,'http://localhost').pathname);let file=path.resolve(root,'.'+pathname);if(!file.startsWith(root+path.sep)&&file!==root){res.writeHead(403);return res.end();}
   try{const stat=await fs.stat(file);if(stat.isDirectory()){if(!pathname.endsWith('/')){res.writeHead(302,{Location:pathname+'/'});return res.end();}file=path.join(file,'index.html');}}catch{}
   const buffer=await fs.readFile(file);res.writeHead(200,{'Content-Type':types[path.extname(file)]||'application/octet-stream','Cache-Control':'no-cache'});res.end(buffer);
  }catch{res.writeHead(404,{'Content-Type':'text/html; charset=utf-8'});res.end(await fs.readFile(path.join(root,'404.html')).catch(()=>Buffer.from('Not found')));}
 });server.listen(port,'0.0.0.0',()=>console.log(`FENG preview: http://localhost:${port}`));return server;
}
if(process.argv[1]===fileURLToPath(import.meta.url))serve();
