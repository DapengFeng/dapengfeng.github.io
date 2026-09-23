(() => {
  const defaults = {
    rust: {compiler:'r1850', name:'rustc 1.85.0', flags:'--edition=2021 --crate-type=bin -O'},
    'c++': {compiler:'g142', name:'GCC 14.2', flags:'-std=c++17 -O2 -Wall -Wextra'},
    python: {compiler:'python312', name:'Python 3.12', flags:'', executor:true}
  };
  const escape = text => String(text).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const pair = (en, zh) => `<span class="i18n"><span lang="en" data-lang="en">${escape(en)}</span><span lang="zh-CN" data-lang="zh">${escape(zh)}</span></span>`;
  const lines = values => Array.isArray(values) ? values.map(v=>String(v.text??'')).join('\n').replace(/\x1b\[[0-9;]*m/g,'') : '';
  function renderOutput(node, text, en, zh) {
    // Godbolt output stays verbatim in every mode. Only our own fallback is localized.
    if (text) node.textContent=text;
    else node.innerHTML=pair(en,zh);
  }
  function sourceText(code) {
    const rows = [...code.querySelectorAll('.code-line')];
    if (!rows.length) return code.textContent;
    return rows.map(row=>{const copy=row.cloneNode(true);copy.querySelectorAll('.line-number').forEach(n=>n.remove());return copy.textContent;}).join('\n');
  }
  document.querySelectorAll('code[data-godbolt]').forEach((code,index) => {
    const language=code.dataset.godbolt, preset=defaults[language];if(!preset)return;
    const compiler=code.dataset.compiler||preset.compiler, name=code.dataset.compilerName||preset.name;
    const flags=code.dataset.compilerOptions??preset.flags;
    const panel=document.createElement('section');panel.className='compiler-check';panel.setAttribute('aria-label','Code editor and results / 代码编辑与结果');
    const editorId=`compiler-editor-${index}`;
    panel.innerHTML=`<label class="compiler-editor-label" for="${editorId}">${pair('Editable code · '+language,'可编辑代码 · '+language)}</label><textarea id="${editorId}" class="compiler-editor" spellcheck="false" autocapitalize="off" autocomplete="off" wrap="off"></textarea><div class="compiler-toolbar"><div><span class="compiler-provider">COMPILER EXPLORER</span><div class="compiler-config">${escape(name)}${flags?` <code>${escape(flags)}</code>`:''}</div></div><div class="compiler-actions"><button type="button" class="compiler-reset">${pair('Reset example','恢复示例')}</button><button type="button" class="compiler-run">${pair('Compile & run','编译并运行')} <span aria-hidden="true">↗</span></button></div></div><p class="compiler-note">${pair('Edit above, then send the current code to Godbolt. Edits last until you reload or switch examples.','在上方修改后，将当前代码发送到 Godbolt。刷新页面或切换示例会恢复示例代码。')}</p><p class="compiler-edited-note compiler-note" hidden>${pair('The article’s explanation describes the original example; results below come from your edited code.','文章说明针对原始示例；下方结果来自你修改后的代码。')}</p><div class="compiler-status" role="status" aria-live="polite">${pair('Not run yet.','尚未运行。')}</div><div class="compiler-results" hidden><section class="compiler-diagnostics"><h4>${pair('Compiler diagnostics','编译诊断')}</h4><pre class="compiler-output" tabindex="0"></pre></section><section class="compiler-runtime" hidden><h4>${pair('Standard output · stdout','标准输出 · stdout')}</h4><pre class="compiler-stdout" tabindex="0"></pre><h4>${pair('Standard error · stderr','标准错误 · stderr')}</h4><pre class="compiler-stderr" tabindex="0"></pre></section><details class="compiler-assembly" hidden><summary>${pair('Generated assembly','生成的汇编')}</summary><pre tabindex="0"></pre></details></div><p class="compiler-note">${pair('A successful run does not prove correctness or memory safety.','运行成功不代表代码始终正确或不存在内存安全问题。')}</p>`;
    const pre=code.closest('pre');pre.after(panel);pre.hidden=true;
    const editor=panel.querySelector('.compiler-editor'),button=panel.querySelector('.compiler-run'),restore=panel.querySelector('.compiler-reset'),status=panel.querySelector('.compiler-status'),results=panel.querySelector('.compiler-results'),output=panel.querySelector('.compiler-output'),assembly=panel.querySelector('.compiler-assembly'),runtime=panel.querySelector('.compiler-runtime');
    editor.dataset.language=language;
    code.compilerEditor=editor;
    let baseline=sourceText(code),source=baseline,controller=null,revision=0;
    editor.value=source;editor.rows=Math.min(24,Math.max(8,source.split('\n').length+1));restore.disabled=true;
    const refreshSyntax=window.CodeSyntax?.mount(editor,language)||(()=>{});
    function change(next){
      refreshSyntax();if(next===source)return;source=next;revision++;controller?.abort();controller=null;button.disabled=false;panel.removeAttribute('aria-busy');panel.removeAttribute('data-state');results.hidden=true;runtime.hidden=true;assembly.hidden=true;
      results.querySelectorAll('pre').forEach(n=>n.textContent='');restore.disabled=source===baseline;panel.querySelector('.compiler-edited-note').hidden=source===baseline;
      status.innerHTML=pair('Code changed. Run this version again.','代码已变更，请重新运行当前版本。');
    }
    function syncExample(){const next=sourceText(code);if(next===baseline)return;baseline=next;editor.value=next;change(next);restore.disabled=true;panel.querySelector('.compiler-edited-note').hidden=true;}
    new MutationObserver(syncExample).observe(code,{subtree:true,childList:true,characterData:true});
    editor.addEventListener('input',()=>change(editor.value));
    restore.addEventListener('click',()=>{editor.value=baseline;change(baseline);editor.focus();});
    button.addEventListener('click',async()=>{
      syncExample();change(editor.value);const current=++revision,snapshot=source;
      controller=new AbortController();const active=controller;let timedOut=false;
      const timer=setTimeout(()=>{timedOut=true;active.abort();},30000);
      button.disabled=true;panel.setAttribute('aria-busy','true');panel.dataset.state='pending';results.hidden=true;runtime.hidden=true;assembly.hidden=true;
      status.innerHTML=pair('Compiling and running on Godbolt…','正在 Godbolt 上编译并运行…');
      try{
        const response=await fetch(`https://godbolt.org/api/compiler/${encodeURIComponent(compiler)}/compile`,{
          method:'POST',credentials:'omit',referrerPolicy:'no-referrer',headers:{'Content-Type':'application/json','Accept':'application/json'},signal:active.signal,
          body:JSON.stringify({source:snapshot,lang:language,allowStoreCodeDebug:false,options:{userArguments:flags,compilerOptions:{skipAsm:false,executorRequest:Boolean(preset.executor)},filters:{execute:true,binary:false,labels:true,directives:true,commentOnly:true,intel:true,demangle:true},executeParameters:{args:[],stdin:''}}})
        });
        if(current!==revision)return;if(!response.ok)throw Error(`HTTP ${response.status}`);
        const data=await response.json();if(current!==revision)return;if(!Number.isInteger(data.code))throw Error('Unexpected API response');
        const execution=data.execResult||(preset.executor?data:null),build=execution?.buildResult||data.buildResult;
        const compiled=(build?.code??(preset.executor?null:data.code))===0;
        const ran=execution?.didExecute===true,exit=execution?.code;
        const passed=compiled&&ran&&exit===0&&!execution.timedOut;
        panel.dataset.state=passed?'passed':'rejected';
        if(!compiled)status.innerHTML=pair(`Compilation failed · exit code ${build?.code??data.code}`,`编译失败 · 退出码 ${build?.code??data.code}`);
        else if(execution?.timedOut)status.innerHTML=pair('Compilation passed · execution timed out','编译通过 · 运行超时');
        else if(ran&&Number.isInteger(exit))status.innerHTML=pair(`Compilation passed · program exit code ${exit}`,`编译通过 · 程序退出码 ${exit}`);
        else status.innerHTML=pair('Compilation passed · program did not run','编译通过 · 程序未运行');
        const diagnostics=[...(preset.executor?[]:[lines(data.stdout),lines(data.stderr)]),lines(build?.stdout),lines(build?.stderr)].filter(Boolean);
        renderOutput(output,[...new Set(diagnostics)].join('\n'),'(No compiler diagnostics)','（无编译诊断）');
        if(execution&&(ran||execution.timedOut||lines(execution.stdout)||lines(execution.stderr))){
          runtime.hidden=false;
          renderOutput(panel.querySelector('.compiler-stdout'),lines(execution.stdout),'(No standard output)','（无标准输出）');
          renderOutput(panel.querySelector('.compiler-stderr'),lines(execution.stderr),'(No standard error)','（无标准错误输出）');
          if(execution.truncated){const notice=document.createElement('span');notice.className='compiler-output-notice';notice.innerHTML=pair('(Output truncated by Godbolt)','（Godbolt 已截断输出）');panel.querySelector('.compiler-stderr').append('\n',notice);}
        }
        const asm=lines(data.asm);if(compiled&&asm){assembly.hidden=false;assembly.querySelector('pre').textContent=asm;}results.hidden=false;
      }catch(error){
        if(current!==revision)return;panel.dataset.state='error';status.innerHTML=timedOut?pair('Godbolt timed out. Retry when ready.','Godbolt 请求超时，可以重试。'):pair('Could not reach Godbolt or read its response. No result; retry.','无法连接 Godbolt 或读取其响应。本次没有结果，请重试。');
        if(error.name==='AbortError')renderOutput(output,'','Request timed out','请求超时');
        else if(error.message==='Unexpected API response')renderOutput(output,'','Unexpected API response','接口返回格式异常');
        else output.textContent=String(error.message);
        results.hidden=false;
      }finally{clearTimeout(timer);if(current===revision){controller=null;button.disabled=false;panel.removeAttribute('aria-busy');}}
    });
  });
})();
